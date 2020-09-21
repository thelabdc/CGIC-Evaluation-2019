## analysis folder:
This folder contains the notebooks and scripts that are used in the analysis. The decision to use both Python and R is largely based on the comfort level of the analysts. The notebooks are listed by the _type_ of analysis, so prepping data for analysis comes first, followed by descriptive analyses, and then by quasi-experimental analysis. For Lab @ DC analysts, the order of these notebooks should actually be run in the order of: 1_, 4_, 5_, 6_, 7_, 8_, 2_, 3_. The reason for this is because the quasi-experimental analyses were our primary task and so this analysis was run first, so files created from running these notebooks are generated out of order. 


- `1_Prep_filter_gun_violent_arrests.ipynb` filters the arrest and crimes data from MPD to generate dataframes for use in the `4_Quasi_Experimental_Analysis_Synthetic_Control.ipynb` notebook.

- `2_Descriptive_Analysis_Clearance_Rates.ipynb` takes a look at the clearance rates for MPD cases with and without NIBIN information. 

- `3_Descriptive_Analysis_Prosecutorial_Outcomes.ipynb` takes a descriptive look at the prosecutorial outcomes for cases which have had a NIBIN lead or NIBIN hit. 

- `4_Quasi_Experimental_Analysis_Synthetic_Control.ipynb` is the main file for running the synthetic control analysis for each of the metrics laid out in the [pre-analysis plan](https://osf.io/q8r5m/). This notebook uses the `synth.py` script from the *cgic_script* folder. Essentially, this file generates rolling averages for every metric for every _PSA within 7D_. Then we generate synthetic rolling averages for every metric, for every day, within the study period, and compare that to the _real_ rolling averages. There are multiple outputs generated from this file: 
    (1) synthetic control plots, which are included in the _outputs_ folder, 
    (2) cleaned files for every metric for running diff-in-diff analysis in R script `7_Quasi_Experimental_Analysis_Diff_in_Diff.R`, and 
    (3) rolling averages for every metric for every PSA for running placebo tests in `5_Data_Prep_Placebo_Tests.ipynb` and `6_Quasi_Experimental_Analysis_Diff_in_Diff.R`. 

- `5_Data_Prep_Placebo_Tests.ipynb` takes the rolling averages generated from `4_Quasi_Experimental_Analysis_Synthetic_Control.ipynb` and prepares the data for `6_Quasi_Experimental_Analysis_Diff_in_Diff.R`. This notebook essentially generates synthetic rolling means for _every PSA in_ _**DC**_, for every metric studied, through multiple loops. This notebook outputs the 192 files to the data/data_placebo folder.

**To run the R scripts that follow:** 

- Open up the R project named `analysis.Rproj` _**first**_. This sets the working directory to be _../CGIC_Evaluation_2019/analysis/_, so you do not need to change the path for where the data lives. 

- `6_Quasi_Experimental_Analysis_Placebo_Tests.R` takes the data generated from `5_Data_Prep_Placebo_Tests.ipynb` to calculate and generate the placebo tests and the error plots that go into our report. These plots are also saved to the _outputs_ folder. 

- `7_Quasi_Experimental_Analysis_Diff_in_Diff.R` takes the cleaned data from `4_Quasi_Experimental_Analysis_Synthetic_Control` and generates a daily count of incidents for each metric. Then we calculate run our difference-in-differences analysis and corresponding confidence intervals, and save them to the _output_ folder. 

- `8_Other_Analysis_Pre_Treatment_Trends.R` generates pre-treatment trend figures to assess the potential vailidity of the parallel trends assumption prior to running our difference-in-differences analysis. The plots are not included in the report itself, but saved to the _output_ folder. 