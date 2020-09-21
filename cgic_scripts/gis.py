"""
Functions for dealing with the GIS data

@author Kevin H. Wilson <kevin.wilson@dc.gov>
"""
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point


"""
EPSG = European Petroleum Survey Group. 
They publish a database of coordinate system information. 
We need to identify a coordinate reference system when reading in geodata
"""

MD_METER_CRS = 2804  # MD State Plane EPSG (This is what DC gov uses)
WSG84_LON_LAT_CRS = 4326  # Lon/Lat EPSG


def read_geofile(filename, xcol, ycol, epsg_crs, **kwargs):
  """
  Read in an Excel or CSV file which contains geodata as two columns
  and convert to a GeoDataFrame.
  
  Args:
    filename (str): The file to read in. Must end in xls, xlsx, or csv
    xcol (str): The column containing the x coordinate
    ycol (str): The column containing the y coordinate
    epsg_crs (int): The EPSG coordinate reference system the data is in
    **kwargs: To be passed to the read_excel/read_csv function of pandas
    
  Returns:
    gpd.GeoDataFrame: The parsed data
  """
  if filename.endswith('xls') or filename.endswith('xlsx'):
    df = pd.read_excel(filename, **kwargs)
  else:
    df = pd.read_csv(filename, **kwargs)
  
  # Filter out bad rows that have non-floats for x/y
  df[xcol] = pd.to_numeric(df[xcol], errors='coerce')
  df[ycol] = pd.to_numeric(df[ycol], errors='coerce')
  df = df[~(df[xcol].isnull() | df[ycol].isnull())]

  # Convert to geodataframe
  gdf = gpd.GeoDataFrame(df, geometry=[Point(lon, lat) for _, (lon, lat) in df[[xcol, ycol]].iterrows()])
  gdf.crs = {'init': 'epsg:{}'.format(epsg_crs)}
  return gdf.to_crs(epsg=MD_METER_CRS)


def sjoin_drop(left, right, how='inner'):
  """
  Do a spatial join between left and right geodataframs and
  drop extraneous columns that geopandas includes.

  Args:
    left (gpd.GeoDataFrame): The left dataframe.
    right (gpd.GeoDataFrame): The right dataframe.
    how (str): How to join the data frames (left, right, inner, outer)

  Returns:
    gpd.GeoDataFrame: gpd.sjoin(left, right) with `index_(left|right)` dropped
  """
  df = gpd.sjoin(left, right, how=how)
  drop_cols = [col for col in ['index_left', 'index_right'] if col in df.columns]
  df.drop(drop_cols, axis=1, inplace=True)
  return df
