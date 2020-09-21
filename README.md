# CGIC-Evaluation-2019
@authors: [Vicky Mei](vicky.mei@dc.gov), [Felix Owusu](felix.owusu@dc.gov), Kevin Wilson

## Introduction
This repo contains scripts and some of the data used to complete the evaluation of the CGIC 2.0 improvements that were piloted in the DC Metropolitan Police Department's (MPD) Seventh District (7D) between the period November 2017 to April 2019. The analysis was pre-registered on the Open Science Framework (OSF). You can read the pre-analysis plan [here](https://osf.io/q8r5m/) and a copy of the final report that was submitted to the Department of Justice (DOJ) Office of Justice Programs (OJP) Bureau of Justice Assistance (BJA) [here](https://osf.io/znwsq/).

#### Requirements:
The primary modules used here (and the one with most issues in terms of versioning) for Python is geopandas, version 0.4.0. We also use numpy, matplotlib, and pandas. The relevant libraries are included at the top of each Python and R scripts and Jupyter notebooks.

## Organization
The repo is organized in four primary folders:
- analysis
- cgic_scripts
- data
- output

## data folder: 
This folder contains all of the data necessary for the CGIC analysis. However, the data available on this repo is strictly publicly available data, i.e. all of the geospatial data that was used in this analysis, and data generated from running our synthetic control analysis, i.e. the real and synthetic rolling averages for every PSA for every metric we are studying. Besides the geojson file for the 2018 Police Service Areas (which we obtained from MPD as it was no longer publicly available), all geospatial files come from ArcGIS's Open Data Hub.  

#### Geospatial data

The links for the open data sets we used are:

* [2019 Police service area polygons](https://opendata.arcgis.com/datasets/db24f3b7de994501aea97ce05a50547e_10.geojson)
* [Advisory Neighborhood Commission polygons](https://opendata.arcgis.com/datasets/fcfbf29074e549d8aff9b9c708179291_1.geojson)
* [Single Member District polygons](https://opendata.arcgis.com/datasets/890415458c4c40c3ada2a3c48e3d9e59_21.geojson)
* [block_groups](https://opendata.arcgis.com/datasets/c143846b7bf4438c954c5bb28e5d1a21_2.geojson)

Here's a quick copy and paste script to grab them:

```bash
wget https://opendata.arcgis.com/datasets/db24f3b7de994501aea97ce05a50547e_10.geojson
wget https://opendata.arcgis.com/datasets/fcfbf29074e549d8aff9b9c708179291_1.geojson
wget https://opendata.arcgis.com/datasets/890415458c4c40c3ada2a3c48e3d9e59_21.geojson
wget https://opendata.arcgis.com/datasets/c143846b7bf4438c954c5bb28e5d1a21_2.geojson
```
MPD updated their police districts as of January 10, 2019. The police service area polygons as linked above are reflected are based on the 2019 update, however, the 2018 police service areas are also included in the folder and used for the analysis. 

#### Datasets that we used to create the report but which are _not_ included in this repository due to the sensitive nature of the data:
- MPD data on crimes occuring in the District between October 2016 to April 2019.
- MPD data on charges for arrests in the District between October 2016 to April 2019.
- Department of Forensic Science (DFS) data on the crimes which have had a NIBIN lead or a hit. 
- Calls for service for sounds of gunshots using data from the Office of Unified Communications (OUC) 
- Data from MPD's gunshot detection technology that uses sensors to detect and alert law enforcement agencies of potential gun incidents in the district (Shotspotter)
- Prosecutorial outcomes data from the U.S. Attorney's Office (USAO)

#### data/data_placebos folder:
This folder contains the output generated from running notebook `4_Quasi_Experimental_Analysis_Synthetic_Control.ipynb`. That is, these are files containing the synthetic and real rolling averages for every metric for every PSA.


## cgic_scripts folder:
This folder contains the helper scripts and functions used in the analysis. Detailed information about each script are as follows:

- `clean_ss_and_cfs.py`: Requires data on CAD events, CAD comments, and ShotSpotter alerts. CAD events, as well as CAD comments, 
generally come in one .csv file per year of events or comments. You should go into the script and edit the filenames and directories (lines 6 through 22) as necessary to navigate to where you have the data saved. Takes calls for service and ShotSpotter data, performs some cleaning functions to fix geolocation data, remove outliers, and keep only CAD events with at least one citizen call associated.
Outputs a cleaned .csv file for calls, and one for alerts, each with time and geolocation data, which is fed into the `4_Quasi_Experimental_Analysis_Synthetic_Control.ipynb` notebook in the _analysis_ folder.

- `gis.py` includes functions for handling with GIS (geopandas) data

- `nnls.py` includes a function for dealing non-negative coefficients for a linear regression 

- `plots.py` includes functions for visualizing DC, synthetic control plots, and bootstrap plots

- `shotspotter.py` includes functions for cleaning Shotspotter data so that events close in location and time are aggregated into a single event.

- `synth.py` includes functions for constructing synthetic controls, including the data prep, t-tests, bootstrapping, and the generation of the synthetic control. This pulls in functions from the `nnls.py` and `plots.py` scripts. We primarily call the `synth.py` script in the `4_Quasi_Experimental_Analysis_Synthetic_Control.ipynb` notebook in the _analysis_ folder.

- `synth_placebo.py` modifies the original synth script to allow for manually including or excluding certain PSAs into the analysis, for the purposes of using the outputs in the placebo tests. We call this script in the `5_Data_Prep_Placebo_Tests.ipynb` notebook.

## analysis folder:
This folder contains the notebooks and scripts that are used in the analysis. The decision to use both Python and R is largely based on the comfort level of the analysts. The notebooks are listed by the _type_ of analysis, so prepping data for analysis comes first, followed by descriptive analyses, and then by quasi-experimental analysis. For Lab @ DC analysts, the order of these notebooks should actually be run in the order of: 1_, 4_, 5_, 6_, 7_, 8_, 2_, 3_. The reason for this is because the quasi-experimental analyses were our primary task and so this analysis was run first, so files created from running these notebooks are generated out of order. 


- `1_Prep_filter_gun_violent_arrests.ipynb` filters the arrest and crimes data from MPD to generate dataframes for use in the `4_Quasi_Experimental_Analysis_Synthetic_Control.ipynb` notebook.

- `2_Descriptive_Analysis_Clearance_Rates.ipynb` takes a look at the clearance rates for MPD cases with and without NIBIN information. 

- `3_Descriptive_Analysis_Prosecutorial_Outcomes.ipynb` takes a descriptive look at the prosecutorial outcomes for cases which have had a NIBIN lead or NIBIN hit. This data cannot be made public, so someone without access to the data cannot run through the code in this file. However, we show the output in these files. 

- `4_Quasi_Experimental_Analysis_Synthetic_Control.ipynb` is the main file for running the synthetic control analysis for each of the metrics laid out in the [pre-analysis plan](https://osf.io/q8r5m/). This notebook uses the `synth.py` script from the *cgic_script* folder. Essentially, this file generates rolling averages for every metric for every _PSA within 7D_. Then we generate synthetic rolling averages for every metric, for every day, within the study period, and compare that to the _real_ rolling averages. There are multiple outputs generated from this file: 
    (1) synthetic control plots, which are included in the _outputs_ folder, 
    (2) cleaned files for every metric for running diff-in-diff analysis in R script `7_Quasi_Experimental_Analysis_Diff_in_Diff.R`, and 
    (3) rolling averages for every metric for every PSA for running placebo tests in `5_Data_Prep_Placebo_Tests.ipynb` and `6_Quasi_Experimental_Analysis_Diff_in_Diff.R`. For every file ending in "_pre", the rows correspond to every day in the pre-intervention period (i.e., the 365 days from 11/1/2016 to 10/31/2017). Similarly, for every file ending in "_post", the rows correspond to every day in the intervenion period (i.e. the 546 days that span between 11/1/2017 and 4/30/2019).      
- `5_Data_Prep_Placebo_Tests.ipynb` takes the rolling averages generated from `4_Quasi_Experimental_Analysis_Synthetic_Control.ipynb` and prepares the data for `6_Quasi_Experimental_Analysis_Diff_in_Diff.R`. This notebook essentially generates synthetic rolling means for _every PSA in_ _**DC**_, for every metric studied, through multiple loops. This notebook outputs the 192 files to the data/data_placebo folder.

**To run the R scripts that follow:** 

- Open up the R project named `analysis.Rproj` _**first**_. This sets the working directory to be _../CGIC_Evaluation_2019/analysis/_, so you do not need to change the path for where the data lives. 

- `6_Quasi_Experimental_Analysis_Placebo_Tests.R` takes the data generated from `5_Data_Prep_Placebo_Tests.ipynb` to calculate and generate the placebo tests and the error plots that go into our report. These plots are also saved to the _outputs_ folder. 

- `7_Quasi_Experimental_Analysis_Diff_in_Diff.R` takes the cleaned data from `4_Quasi_Experimental_Analysis_Synthetic_Control` and generates a daily count of incidents for each metric. Then we calculate run our difference-in-differences analysis and corresponding confidence intervals, and save them to the _output_ folder. 

- `8_Other_Analysis_Pre_Treatment_Trends.R` generates pre-treatment trend figures to assess the potential vailidity of the parallel trends assumption prior to running our difference-in-differences analysis. The plots are not included in the report itself, but saved to the _output_ folder. 