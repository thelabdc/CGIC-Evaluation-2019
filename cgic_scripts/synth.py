"""
Some functions for constructing synthetic controls

@author Kevin H. Wilson <kevin.wilson@dc.gov>
@author Vicky Mei <vicky.mei@dc.gov>
"""
import numpy as np
import pandas as pd
from scipy import stats as st
from sklearn import linear_model

from . import nnls, plots


def aggregate_by_col(df, col, num_days=30,
                     min_date='2016-10-01', max_date='2019-05-01',
                     time_col='event_time', measure = 'RM'):
  """
  Return a data frame which contains the rolling mean of the number of events per day
  in `num_days` windows starting at `min_date` grouped by `col`

  Args:
    df (pd.DataFrame): This is the data frame to aggregate
    col (str): The name of the column we're aggregating by
    num_days (int): The window we're averaging over
    min_date (str): The minimum date for aggregation in YYYY-MM-DD format
    max_date (str): The maximum date for aggregation in YYYY-MM-DD format
    time_col (str): The name of the column representing time
    measure (str): Default is 'RM' for 'rolling mean'. Else, input 'M' or 'W' to get the 
                   monthly or weekly average

  Returns:
    pd.Series: A series with a two-level index. The first level is `col`, the
      second level is the date, and the value is rolling mean
  """
  df = df[(df[time_col] >= min_date) & (df[time_col] < max_date)]
  count_df = df.groupby([col, pd.Grouper(key=time_col, freq='D')]).size().unstack(col)
  idx = pd.date_range(min_date, max_date, closed='left')
  idx.name = 'event_time'
  count_df = count_df.reindex(idx).fillna(0)
  if measure == 'RM':
    count_df = count_df.rolling(window=30).mean()
    return count_df.loc[pd.datetime(2016, 11, 1):].stack(col).reorder_levels([1, 0]).sort_index()
  else: 
    return count_df.groupby(pd.Grouper(freq=measure)).mean().stack(col).reorder_levels([1, 0]).sort_index().dropna()


def _prep_data(df, col, min_time, max_time, control_min_val, control_max_val):
  """
    Return two arrays with the rolling means of the number of events per day 
    for the control group (between min_ and max_ times and between control_min_val and
    control_max_val) and the treatment group (greater than max_time and control_max_val)
    
    Args: 
        df (pd.DataFrame): The dataframe containing the rolling means by geography (eg PSA) and day
                           (pass through aggregate_by_col first)
        col (str): The name of the region we want to create a synthetic control for
        min_time (str): The year before the intervention started in YYYY-MM-DD format
        max_time (str): The date the intervention started in YYYY-MM-DD format
        control_min_val (int): Control group PSA starting val (eg 101 for psa)
        control_max_val (int): Control group PSA ending value (eg 700 for psa)
        
    """
  
  on_time = df[(df.event_time >= min_time) & (df.event_time < max_time)]
  control_rows = (on_time[col] >= control_min_val) & (on_time[col] < control_max_val)
  treatment_rows = on_time[col] >= control_max_val

  A = on_time[control_rows].pivot_table('rm', col, 'event_time').sort_index().fillna(0).values
  b = on_time[treatment_rows].pivot_table('rm', col, 'event_time').sort_index().fillna(0).values

  return A, b



def perform_control(ddf, col='event_psa', p=0.75,
                    fig=None, ax=None, plot_all_periods = True,
                    soft_convex=False, random_state=None,
                    control_min_val=101, control_max_val=700,
                    control_min_time='2016-11-01', control_max_time='2017-11-01',
                    treatment_min_time='2017-11-01', treatment_max_time='2019-05-01',
                    bootstrap = False, path_title = ''):
  """
  What does this do?
  
  Args:
      ddf (pd.DataFrame): The dataframe containing the rolling means by geography (eg PSA) and day
                           (pass through aggregate_by_col first)
      col (str): The name of the region we want to create a synthetic control for
      p (float): The proportion of control rows to use in construction of synthetic control
      fig (plt.figsize|None): The size of the figure to construct. Ignored if `ax` is not None
      ax (plt.Axes|None): The axes to plot on or None if we should create new axes
      soft_convex (bool): True if we want to allow for extrapolation
      random_state (int|None): True if we want to set a seed to replicate results,
      control_min_val (int): The starting PSA to use as our control group
      control_max_val (int): The starting PSA of our treatment group,
      control_min_time (str): The start date of the control period in 'YYYY-MM-DD' format
      control_max_time (str): The end date of the control period in 'YYYY-MM-DD' format
      treatment_min_time (str): The start date of the treatment period in 'YYYY-MM-DD' format 
      treatment_max_time (str): The end date of the treatment period in 'YYYY-MM-DD' format
      multiple_tests (bool): Determine if you want to conduct t-tests for every PSA in 7D. Default is False
      
  Returns:
      result: Tuple containing the t-statistic and the p-value result from a two-sided t-test
      residual: Array containing the the differences between the observed and predicted values of the data for 7D.  
      len_test: Int containing the length of the differences in rolling means between the synthetic and real 7D
      weights: N-Dimensional Array (control x treat) containing the weights contributed by every control group PSA to the 
               construction of the synthetic 7D
  
  """    
    
    
  random_state = random_state or np.random

  # Prep data for the control period. 
  A1, b1 = _prep_data(ddf, col, control_min_time, control_max_time, control_min_val, control_max_val)
 
  # If p < 1, then choose a random collection of rows to keep for the training
  # This allows us to test robustness
  # Choose a random sample from the number of control group geographies (eg PSAs) available based on p percentage
  if p < 0.99:
    keep_rows = random_state.choice(A1.shape[0], size=int(A1.shape[0] * p), replace=False)
    A1 = A1[keep_rows, :]

  # Add the relaxed convex condition
  # np.hstack() takes a sequence of arrays and stacks them horizontally to make a single array.  
  if soft_convex:
    A1 = np.hstack([A1, np.ones((A1.shape[0], 1))])
    b1 = np.hstack([b1, np.ones((b1.shape[0], 1))])

  # Optimize with a non-negative least squares
  weights = [] 
  residual = []
  models = []
  # for every police district in the treated geography (PSA), fit a non-negative least squares model
  for i in range(b1.shape[0]):
    model = nnls.NNLS(normalize=False, fit_intercept=False)
    model.fit(X = A1.T, y = b1[i, :]) #Fit using control as training data (A1), with our target values (b1)

    #Append results
    weights.append(model.coef_[:, np.newaxis])
    residual.append(model.residues_)
    models.append(model)

  # np.hstack() takes a sequence of arrays and stacks them horizontally to make a single array.  
  weights = np.hstack(weights)

  # Prep data for the treatment period
  A2, b2 = _prep_data(ddf, col, treatment_min_time, treatment_max_time,
                    control_min_val, control_max_val)
  if p < 0.99:
    A2 = A2[keep_rows, :]
       
  # if fig, plot the data
  if fig:
    if plot_all_periods:
        # see function in other py file
        plots.synth_plot_all_periods(A1, b1, A2, b2, weights, fig=fig, soft_convex=soft_convex, path_title=path_title) 
    else:
        plots.synth_plot(A2, b2, weights, fig=fig, soft_convex=soft_convex)
        

  # Synthetic Control Results for PRE-Intervention Period   
  synthetic_control_pre = A1.T @ weights # Apply the weights to our new data in the treatment period
  
  # Get the difference between the real and the synthetic 7D
  test_pre = (b1.T-synthetic_control_pre)[::30, :]
    
  # ---------------------------Synthetic Control Results for Intervention Period ---------

  # synthetic_control = np.hstack([model.predict(A.T)[:, np.newaxis] for model in models])
  synthetic_control = A2.T @ weights # Apply the weights to our new data in the treatment period
  
  # Get the difference between the real and the synthetic 7D
  test = (b2.T - synthetic_control)[::30, :]

  #Initialize to append our results
  result_pre = []
  result = []
  
  if test.shape[1] > 1:  #if we're doing PSA level
      for i in range(test.shape[1]):
        #T-test for pre-intervention period
        t_test_pre = st.ttest_1samp(test_pre[:,i], 0)
        result_pre.append(t_test_pre)
        #T-test for intervention period
        t_test_post = st.ttest_1samp(test[:,i], 0)
        result.append(t_test_post)
  else: #if we're doing district level
      test_overall = test.flatten()
      result_overall = st.ttest_1samp(test_overall, 0)
      pre_period_result = st.ttest_1samp(test_pre.flatten(), 0)
      print('Results for Pre-Intervention Period:', pre_period_result)
      print('Result after intervention:', result_overall)

  if bootstrap == False:
    if test.shape[1] > 1:
        print('Pre')
        for i in range(test.shape[1]):
            print('70{}'.format(i+1), result_pre[i])

        print('Post') 
        for i in range(test.shape[1]):
            print('70{}'.format(i+1), result[i])

        test_overall = test.flatten()
        result_overall = st.ttest_1samp(test_overall, 0)
        print('\n')
        print('Overall results:', result_overall)
    else:
        pass
              

  if bootstrap:
    return result
  else:
    return result_pre, result, residual, weights, synthetic_control, b2.T, synthetic_control_pre, b1.T


def bootstrap_control(ddf, p=0.75, random_state=None, num_to_bootstrap=100, fig=None, bootstrap = True):
  """
  Bootstrap many different synthetic control experiments from the passed data frame.
  This is done by taking the proportion `p` of the *control* rows in `ddf` and
  forming a synthetic control for all the treatment rows.
  """
  bootstrapped = [perform_control(ddf, p=p, random_state=random_state, bootstrap = True)
                  for _ in range(num_to_bootstrap)]
  return bootstrapped
