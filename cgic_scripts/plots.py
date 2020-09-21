"""
Special plots that make it easy to visualize DC

@author Kevin H. Wilson <kevin.wilson@dc.gov>
@author Vicky Mei <vicky.mei@dc.gov>
"""
import itertools as its #For efficient looping

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

#The following two are for label warnings
import sys
import warnings

#colors for plotting
lab_blue = '#2b4888'
lab_pink = '#de4057'
lab_grey = '#595959'

def outline_plot(data_df, data_col, outline_df, ax=None, alpha=0.2, figsize=(10, 10)):
  """
  Plot `data_col` of `data_df` on top of the polygons described by `outline_df`.
  Perform this on `ax` unless ax in None, in which case, create a new figure
  and plot there.

  Args:
    data_df (gpd.GeoDataFrame): The point data to plot
    data_col (str): Which col of `data_df` to plot
    outline_df (gpd.GeoDataFrame): The outline data to plot
    ax (plt.Axes|None): The axes to plot on or None if we should create new axes
    alpha (float): The transparency for plotting points from `data_df`
    figsize (tuple[int, int]): The size of the figure to construct. Ignored if `ax` is not None

  Returns:
    (plt.Figure, plt.Axes)|None: If the ax was None, return the figure and axes we created
  """
  was_ax_null = not bool(ax)
  if not ax:
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(1, 1, 1)
  
  outline_df.plot(ax=ax, color='white', edgecolor='black')
  data_df.plot(ax=ax, column=data_col, alpha=alpha)

  if was_ax_null:
    return fig, ax


def bootstrap_plot(results, fig=None):
    
    """
    Plot the bootstrapped distribution of pvalues and the bootstrapped distribution of the t-statistic 

    Args:
    results (nd.array): The results from running a t-test using scipy.stats.ttest_1samp()  
                         containing the t-statistic and the p-value for every PSA for every bootstrap 
    fig (plt.Axes|None): The axes to plot on or None if we should create new axes
        
    Returns:
    (plt.Figure, plt.Axes)|None: If the ax was None, return the figure and axes we created
    """  
    #Convert to a dataframe
    results = pd.DataFrame(results)
    fig1, axes = plt.subplots(nrows=2, ncols=4, figsize = (20,8), sharey=True)
    fig1.subplots_adjust(hspace=.3) #give more space between subplots
    fig1.suptitle('"Bootstrapped" distribution of t statistic')

    #start the iterator
    iterate = 0
    for i in range(2):
        for j in range(4):
            ax = axes[i,j]
            tstat, _ = list(zip(*results.iloc[:,iterate]))
            ax.hist(tstat, bins=20, density=True, color = lab_grey)
            ax.set_title('PSA {}'.format(700 + iterate + 1))
            ax.set_xlabel('t-statistic')
            ax.set_yticks([])
            iterate += 1
    print()        

    fig2, axes = plt.subplots(nrows=2, ncols=4, figsize = (20,8), sharey=True)
    fig2.subplots_adjust(hspace=.3) #give more space between subplots
    fig2.suptitle('"Bootstrapped" distribution of p-values')        
    iterate = 0
    for i in range(2):
        for j in range(4):
            ax = axes[i,j]
            _, pval = list(zip(*results.iloc[:,iterate]))
            ax.hist(pval, np.arange(0, 1.01, 0.05), density=True, color=lab_blue)
            ax.set_title('PSA {}'.format(700 + iterate + 1))
            ax.set_xlabel('p-values')
            ax.set_yticks([])

            iterate += 1

           
def bootstrap_plot_district(results, fig=None):
    
    """
    Plot the bootstrapped distribution of pvalues and the bootstrapped distribution of the t-statistic 

    Args:
    results (nd.array): The results from running a t-test using scipy.stats.ttest_1samp()  
                         containing the t-statistic and the p-value for every PSA for every bootstrap 
    fig (plt.Axes|None): The axes to plot on or None if we should create new axes
        
    Returns:
    (plt.Figure, plt.Axes)|None: If the ax was None, return the figure and axes we created
    """  
    fig = fig or plt.figure(figsize=(10, 5))
  
    ax1, ax2 = fig.subplots(1, 2)

    ax1.hist([b[0][1] for b in results], np.arange(0, 1.01, 0.05), density=True, color = lab_blue)
    ax1.set_title('"Bootstrapped" distribution of p-values')
    ax1.set_xlabel('p-value')
    ax1.set_yticks([])
  
    ax2.hist([b[0][0] for b in results], bins=20, density=True, color = lab_grey)
    ax2.set_title('"Bootstrapped" distribution of t-statistic')
    ax2.set_xlabel('t-statistic')
    ax2.set_yticks([])

    return fig

def synth_plot(A, b, weights, fig=None, soft_convex=False):
    
  """
  Plot the rolling means for the PSAs in 7D for the synthetic 7D against the real 7D. 

  Args:
    A (np.ndarray): N-dimensional Array of rolling means for every PSA in the 
                    control group (1D-6D). Typically 40 by X number of days 
    b (np.ndarray): N-dimensional Array of rolling means for every PSA in the 
                    treatment group (7D). Typically 8 by X number of days
    weights (np.ndarray): N-dimensional Array of weights contributed by the PSAs in the 
                          control group. Typically 40 (control group PSAs) x 8 (7D PSAs). 
    fig (plt.Axes|None): The axes to plot on or None if we should create new axes
    
  Returns:
    (plt.Figure, plt.Axes)|None: If the ax was None, return the figure and axes we created
  """  
  # Suppress label warnings
  if not sys.warnoptions:
    warnings.simplefilter("ignore")

  # Colors to plot
  lab_blue = '#2b4888'
  lab_pink = '#de4057'
    
  fig = fig or plt.figure(figsize=(16, 8))

  # The synthetic control as generated from our model
  # Array of Weights for X number of Days by 8 PSAs (in 7D)  
  synthetic_control = A.T @ weights
    
  if soft_convex:
    to_plot_x = synthetic_control[:-1, :]
    to_plot_y = b[:, -1]
  else:
    to_plot_x = synthetic_control #our synthetic 7D
    to_plot_y = b # our real 7D

  # If there is only 1 synthetic control area, eg District instead of PSA, 
  # only do a 1 plot subplot. Otherwise 2x4 for the 8 PSAs in 7D
  if to_plot_y.shape[0] == 1:
    subplots = [fig.add_subplot(1, 1, 1)]
  else:
    subplots = its.chain(*fig.subplots(2, 4, sharex=True, sharey=True))
    
  for i, ax in enumerate(subplots):
    real = ax.plot(to_plot_y[i, :], c=lab_blue) # Real 7D 
    synthetic = ax.plot(to_plot_x[:, i], c=lab_pink, linestyle='dashed') # Synthetic 7D
    ax.set_xticklabels([])
    ax.set_xticks([])
    if to_plot_y.shape[0] == 1:
        ax.set_title('7D'.format(700 + i + 1))
    else:
        ax.set_title('PSA {}'.format(700 + i + 1))
        

  fig.text(0.08, 0.5, 'Rolling 30-day rate', va='center', rotation = 'vertical')
  fig.legend([real, synthetic], labels = ['Real', 'Synthetic'], loc = 'upper right', frameon = False, borderaxespad=1)

    
def synth_plot_all_periods(A1, b1, A2, b2, weights, fig=None, soft_convex=False, path_title = ''):
    
  """
  Plot the rolling means for the PSAs in 7D for the synthetic 7D against the real 7D. 

  Args:
    A (np.ndarray): N-dimensional Array of rolling means for every PSA in the 
                    control group (1D-6D). Typically 40 by X number of days 
    b (np.ndarray): N-dimensional Array of rolling means for every PSA in the 
                    treatment group (7D). Typically 8 by X number of days
    weights (np.ndarray): N-dimensional Array of weights contributed by the PSAs in the 
                          control group. Typically 40 (control group PSAs) x 8 (7D PSAs). 
    fig (plt.Axes|None): The axes to plot on or None if we should create new axes
    soft_convex (bool): Allow for extrapolation? True/False
    
  Returns:
    (plt.Figure, plt.Axes)|None: If the ax was None, return the figure and axes we created
  """  
  # Suppress label warnings
  if not sys.warnoptions:
    warnings.simplefilter("ignore")

  # Colors to plot
  lab_blue = '#2b4888'
  lab_pink = '#de4057'
  lab_grey = '#595959'
    
  fig = fig or plt.figure(figsize=(16, 8))

  # The synthetic control as generated from our model
  # Array of Weights for X number of Days by 8 PSAs (in 7D)
  # _pre is the pre-intervention period; _post is the post-intervention period
  synthetic_pre = A1.T @ weights
  synthetic_post = A2.T @ weights
  
  synthetic_control = np.append(synthetic_pre, synthetic_post, axis = 0)
    
  #Appending our real 7Ds together
  real_7D = np.append(b1.T, b2.T, axis = 0)
  
  if soft_convex:
    to_plot_x = synthetic_control[:-1, :]
    to_plot_y = b[:, -1]
  else:
    to_plot_x = synthetic_control #our synthetic 7D
    to_plot_y = real_7D # our real 7D

  # If there is only 1 synthetic control area, eg District instead of PSA, 
  # only do a 1 plot subplot. Otherwise 2x4 for the 8 PSAs in 7D
  if to_plot_y.shape[1] == 1:
    subplots = [fig.add_subplot(1, 1, 1)]
  else:
    subplots = its.chain(*fig.subplots(4, 2, sharex=True, sharey=True))
    
  for i, ax in enumerate(subplots):
    real = ax.plot(to_plot_y[:, i], c=lab_blue) # Real 7D 
    synthetic = ax.plot(to_plot_x[:, i], c=lab_pink, linestyle='dashed') # Synthetic 7D
    ax.axvline(x = len(synthetic_pre)-1, color = lab_grey, linestyle = '--', linewidth = .9) #end of control pd
    ax.set_xticklabels([])
    ax.set_xticks([])
    if to_plot_y.shape[1] == 1:
        ax.set_title('7D')
    else:
        ax.set_title('PSA {}'.format(700 + i + 1))
        
  fig.text(0.08, 0.5, 'Rolling 30-day rate', va='center', rotation = 'vertical')
  fig.legend([real, synthetic], labels = ['Real', 'Synthetic'], frameon = False, ncol = 3, fontsize = 'medium')
  fig.subplots_adjust(top = .94)

  fig.savefig(path_title+'.pdf')

