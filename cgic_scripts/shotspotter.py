"""
Clean duplicate ShotSpotter events.

@author Kevin H. Wilson <kevin.wilson@dc.gov>
"""
from collections import Counter
from datetime import timedelta

import geopandas as gpd
import networkx as nx
import pandas as pd
from shapely.geometry import Point

from . import gis #See .py file in cgic folder


def clean_duplicates(df, spatial_buffer=50, temporal_buffer=2, verbose=False):
  """
  The ShotSpotter data seems to have a lot of duplicate *events*, i.e.,
  events that are very close in space and time to each other. This function
  takes a data frame filled with shot spotter data, and aggregates the events
  into those falling within `spatial_buffer` meters of each other and `temporal_buffer`
  minutes.

  The returned data frame will have the following columns as described:
    * `event_id`: If the event represents a grouped event, this will be negative.
      Else it will be the original event id.
    * `event_time`: The *minimum* time of all events in the group
    * `mean_rounds`: The mean number of rounds of each event in the group
    * `total_rounds`: The total number of rounds in the group
    * `geometry`: Will be the centroid of all the event points in the group.
      Note that this is *not* weighted by the number of rounds or any other
      factor.

  Args:
    df (pd.DataFrame): A data frame consisting of ShotSpotter data. It must have
      the following columns: `event_time`, `event_id`, `geometry`, `rounds`.
    spatial_buffer (float): The number of units to buffer each point in `df` when
      trying to find duplicates. This should be in the same units as the `geometry`
      column of `df`.
    temporal_buffer (float): The number of *minutes* in which to consider incidents
      the same.
    verbose (bool): If true, print statistics about the deduplication during the
      deduplication.

  Returns:
    pd.DataFrame: A data frame with the columns as described above
  """
  if verbose:
    print('Before deduplicating, there were {} events'.format(len(df)))

  ## Determine events that are close in space/location

  # First, draw a polygon around every point for every event in the dataframe.
  # buffered will be a gdf of event ids, event times, and polygons
  buffered = gpd.GeoDataFrame({
    'right_event_id': df.event_id.values,
    'right_event_time': df.event_time.values
  }, geometry=df.buffer(spatial_buffer).values)

  # spatial join the shotspotter data (left) to the polygon data in the buffered gdf (right)
  # works if the point in the shotspotter data is located within the polygon in the buffered gdf
  close_events = gis.sjoin_drop(df[['event_id', 'event_time', 'geometry']], buffered)
    
  # bools for whether the shotspotter data id matches the buffered gdf id  
  close_filter_space = close_events.event_id != close_events.right_event_id

  # Determine events that are close in time
  between_time = timedelta(minutes=temporal_buffer)
  # True if the difference between the event times are within 2 minutes of each other
  close_filter_time = ((close_events.event_time - close_events.right_event_time < between_time) &
                       (close_events.right_event_time - close_events.event_time < between_time))

  # Returns a series of booleans if both close in location and in time. 
  close_filter = close_filter_space & close_filter_time

  # Restrict to events that are close in both space and time
  only_close_events = close_events[close_filter]

  if verbose:
    print('We found that {} events were close in space and time to another event'.format(
      only_close_events.event_id.nunique()))

  # Find all connected components of the closeness graph
  # Will return the events that are connected, and which events are connected

  # create an empty graph with no nodes and no edges
  graph = nx.Graph() 
  # Add unique event_ids as our nodes
  graph.add_nodes_from(only_close_events.event_id.unique())
  # Add edges to our nodes the event ids
  for _, row in only_close_events.iterrows():
    graph.add_edge(row.event_id, row.right_event_id)

     
  if verbose:
    sizes = Counter(len(c) for c in nx.connected_components(graph))
    print('The number of connected components of given size:')
    for size, count in sorted(sizes.items(), key=lambda x: x[0]):
      print('\t{}\t{}'.format(size, count))

  # Annotate each event with a group id
  # Multiple events will have the same group id if they are close in time/space
  group_ids = [
    (-group_id, event_id) for group_id, component in enumerate(nx.connected_components(graph), 1)
    for event_id in component
  ]
  group_df = pd.DataFrame.from_records(group_ids, columns=['group_id', 'event_id'])
  with_group_df = pd.merge(df, group_df, how='left', on='event_id')
  with_group_df['group_id'] = with_group_df.group_id.fillna(with_group_df.event_id)

  # For all events that are close in time & space, set the event time to be the first incident, 
  # and the location and rounds to be the average location/rounds of the two events
  to_group = with_group_df[with_group_df.group_id < 0]
  min_time = to_group.groupby('group_id').event_time.min()
  mean_location = to_group.groupby('group_id').geometry.agg(lambda x: Point(x.x.mean(), x.y.mean()))
  mean_rounds = to_group.groupby('group_id').rounds.mean()
  total_rounds = to_group.groupby('group_id').rounds.sum()

  grouped_events_df = min_time.to_frame(name='event_time')\
    .join(mean_rounds.to_frame(name='mean_rounds'))\
    .join(total_rounds.to_frame(name='total_rounds'))\
    .join(mean_location.to_frame(name='geometry'))
  grouped_events_df = gpd.GeoDataFrame(grouped_events_df)

  grouped_events_df.crs = {'init': 'epsg:{}'.format(gis.MD_METER_CRS)}
  grouped_events_df.index.name = 'event_id'
  grouped_events_df = grouped_events_df.reset_index()

  # For events that are not connected to another, add the relevant columns for easy concatentation
  single_events_df = with_group_df[with_group_df.group_id >= 0].copy()
  single_events_df['total_rounds'] = single_events_df.rounds
  single_events_df['mean_rounds'] = single_events_df.rounds

  if verbose:
    print('Retained events: {}'.format(len(grouped_events_df) + len(single_events_df)))

  return pd.concat([
    single_events_df[['event_id', 'event_time', 'mean_rounds', 'total_rounds', 'geometry']],
    grouped_events_df[['event_id', 'event_time', 'mean_rounds', 'total_rounds', 'geometry']],
  ])
