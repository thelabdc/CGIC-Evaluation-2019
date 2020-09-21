"""
This modifies the code from synth.py to allow for manual 
inclusion and exclusion of PSAs from the construction of 
the synthetic control

"""
import numpy as np
import pandas as pd
from scipy import stats as st
from sklearn import linear_model

from . import nnls, plots


def _prep_data(df, col, min_time, max_time, treatment_psas, exclude_psas):
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
        treatment_psas (list): List of ints representing PSAs we want to conduct a synthetic control for
        exclude_psas (list): List of ints representing PSAs we want to exclude from the control group
        
    """
  
  on_time = df[(df.event_time >= min_time) & (df.event_time < max_time)]

  # Find the indices of the PSAs in the treatment group
  treatment_rows = on_time[col].isin(treatment_psas)
  # Find the indices of PSAs in the control group, i.e. NOT in the treatment, and excluding 
  control_rows = ~on_time['event_psa'].isin(treatment_psas) & ~on_time['event_psa'].isin(exclude_psas) 
          
  #A is the control, aka the rows NOT in the treatment rows; b = treatment rows
  A = on_time[control_rows].pivot_table('rm', col, 'event_time').sort_index().fillna(0).values
  b = on_time[treatment_rows].pivot_table('rm', col, 'event_time').sort_index().fillna(0).values

  return A, b



def perform_control(ddf, 
                    col='event_psa', p=0.75,
                    fig=None, ax=None, plot_all_periods = True,
                    soft_convex=False, random_state=None,
                    treatment_psas=[], exclude_psas = [],
                    control_min_time='2016-11-01', control_max_time='2017-11-01',
                    treatment_min_time='2017-11-01', treatment_max_time='2019-05-01',
                    bootstrap = False):
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
      treatment_psas (list): A list of PSAs (ints) we want to construct a synthetic control for
      exclude_psas (list): A list of PSAs (ints) we want to exclude from the control group
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
  A1, b1 = _prep_data(ddf, col, control_min_time, control_max_time, treatment_psas, exclude_psas)
 
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
  A2, b2 = _prep_data(ddf, col, treatment_min_time, treatment_max_time, treatment_psas, exclude_psas)
  if p < 0.99:
    A2 = A2[keep_rows, :]
       
  # if fig, plot the data
  if fig:
        plots.synth_plot(A2, b2, weights, fig=fig, soft_convex=soft_convex)
        

  # Synthetic Control Results for PRE-Intervention Period   
  synthetic_control_pre = A1.T @ weights # Apply the weights to our new data in the treatment period
  
  # Get the difference between the real and the synthetic 7D
  test_pre = (b1.T-synthetic_control_pre)[::30, :]
    
  # ---------------------------Synthetic Control Results for Intervention Period ---------

  # synthetic_control = np.hstack([model.predict(A.T)[:, np.newaxis] for model in models])
  synthetic_control = A2.T @ weights # Apply the weights to our new data in the treatment period
  
  # Get the difference between the real and for every 30 step - to feed into ttest
  test = (b2.T - synthetic_control)[::30, :]
    
  #Get the difference in rolling means between synth and real (for plotting)
  pre_plot = (synthetic_control_pre - b1.T).T
  post_plot = (synthetic_control - b2.T ).T  
  control_less_real = np.concatenate((pre_plot, post_plot), axis = 1) #control (predicted) less real pre + post period
  
  #Initialize to append our results
  result_pre = []
  result = []
  
  for i in range(test.shape[1]): #t-test by PSAs.
    #T-test for pre-intervention period
    t_test_pre = st.ttest_1samp(test_pre[:,i], 0)
    result_pre.append(t_test_pre)
    #T-test for intervention period
    t_test_post = st.ttest_1samp(test[:,i], 0)
    result.append(t_test_post)
  
  if bootstrap == False:
    test_overall = test.flatten()
    result_overall = st.ttest_1samp(test_overall, 0)
              
  if bootstrap:
    return result
  else:
    return synthetic_control, b2.T, synthetic_control_pre, b1.T, control_less_real  
 