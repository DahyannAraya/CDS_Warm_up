![](Logos.jpg)



Welcome to our UNSEQUA documentation, integrating robust impact assessment tools with the valuable climate data from the Copernicus Climate Change Service (C3S) with CLIMADA ETH. We provide two comprehensive notebooks that serve as a foundational gateway into the world of climate impact analysis.

## Installation and Data Access Documentation - Integration of CLIMADA Impact Assessment with Copernicus Climate Change Service (C3S) data

**This first book provides guidance on installing CLIMADA and the Copernicus API for accessing data from the Copernicus Data Store (CDS). It also demonstrates how to explore, select, and download data from the CDS:**

- **How to Install CLIMADA the Copernicus API**:Instructions on how to proceed with the installation. 
- **How to Search for Data**: Guidelines on locating the required data sets.
- **How to Access Data**: Steps to gain access to data once found.
- **Using Copernicus API with CLIMADA ETH**: Instructions on utilizing the Copernicus API within the context of CLIMADA ETH for environmental data analysis.


## Warm-Up documentation - Integration of CLIMADA Impact Assessment with Copernicus Climate Change Service (C3S) data

**This initial exercise is designed as a framework to introduce users to Copernicus data and CLIMADA. Through this exercise, users will learn:**

- **Introduction to CLIMADA ETH for Impact Analysis**: An overview of how CLIMADA ETH can be used for assessing environmental impacts.
- **Data Manipulation**: How to process data to transform it into a hazard object suitable for analysis.
- **Exploring the Hazard Object**: Detailed guide on analyzing the hazard object to extract meaningful insights.
- **Accessing Ready-to-Use Exposure Data (LitPop)**: Procedure to access LitPop, a ready-to-use exposure data set provided by CLIMADA.
- **Accessing Custom Exposure Data**: Instructions on how to search for and use your own exposure data for analysis.
- **Developing an Impact Function**: Steps to create a function to calculate the potential impact based on the data.
- **Calculating the Impact**: Methodology for computing the impact using the developed function.
- **Exploring Impact Results**: Techniques for analyzing and understanding the results of the impact assessment.
- **Saving Analysis Results**: Guidance on how to save the results of the analysis in both tabular and spatial formats for further use.


## Scenario documentation - Integration of CLIMADA Impact Assessment with Copernicus Climate Change Service (C3S) data
**This scenario exercise introduced users to analyzing heatwave impacts under future climate scenarios using CLIMADA and data from the Copernicus Climate Data Store (CDS). Through this exercise, users learned:**

- **Exploring Climate Scenarios**: How to work with climate change scenarios like RCP 8.5 to understand their effects on heatwave trends and impacts on populations.
- **Accessing and Downloading Data**: Steps to search for and download scenario-based data from the Copernicus Climate Data Store, focusing on variables related to heatwaves and population exposure.
- **Setting up the Heatwave Hazard**: How to define and set up the heatwave hazard by integrating climate scenario data, enabling long-term projections of heatwave frequency and intensity.
- **Estimating Population Exposure**: A guide to estimating the population exposed to future heatwaves using socio-economic pathways such as SSP5, highlighting how demographic changes influence risk.
- **Defining Vulnerability and Impact Functions**: Detailed instructions on creating custom impact functions to calculate the potential impacts of heatwaves, considering both climate and population data.
- **Calculating Impacts**: Methodology for calculating the impacts of heatwaves across different decades, assessing how future climate conditions will affect vulnerable populations.
- **Visualizing and Analyzing Results**: Techniques for visualizing the progression of heatwave impacts over time, from 2010 to 2080, allowing users to explore how risks intensify under future climate scenarios.
- **Saving Results for Further Use**: Guidance on saving analysis results in both NetCDF and GeoTIFF formats for use in GIS tools, supporting further spatial analysis and adaptation planning.


## Uncertainty and Sensitivity documentation - Integration of CLIMADA Impact Assessment with Copernicus Climate Change Service (C3S) data

- **Data Extraction**: Retrieved and processed the European heatwave dataset from Copernicus for climate scenarios RCP 8.5.
- **Heatwave Visualization**: Illustrated annual heatwave days for 2010 and 2080, highlighting expected increases in frequency.
- **Uncertainty Analysis**: Employed CLIMADA's unsequa module to assess the impact of uncertain parameters on model outcomes.
- **Impact Function Analysis**: Evaluated the vulnerability function for heatwaves, using a sigmoid function for mean damage degree (MDD) calculation.
- **Sampling Techniques Exploration**: Investigated various SALib sampling methods in CLIMADA for sensitivity analysis.
- **Uncertainty Visualization**: Displayed uncertainty in heatwave impacts on populations through distribution plots.
- **Impact Values Comparison**: Contrasted deterministic impacts with uncertain outcomes to underscore variability in annual impacts and frequency curves.


## Functions

Within the context of this project, specialized functions have been developed to assist users in exporting their findings, particularly the Annual Expected Impact (AEI). AEI refers to the anticipated impact in an average year for each exposure point, available in NetCDF or GeoTIFF formats. This feature is particularly beneficial for integrating CLIMADA's outcomes with other geospatial platforms or software for extended analysis. Moreover, these functionalities facilitate quicker and more straightforward visual communication of the results.

## Data Description

### Heatwave Hazard Analysis

**Objective**: In this first exercise, our goal is to analyze the heatwave hazard in Europe using the heatwave days dataset from Copernicus. This dataset is instrumental in understanding the potential heatwave impacts on the region, informed by climate change projections.

**Dataset Overview**:
- **Source**: Copernicus
- **Climate Scenario**: RCP 8.5 representing different greenhouse gas concentration trajectories.
- **Data Type**: Bias-adjusted EURO-CORDEX dataset.
- **Temporal Coverage**: 1971 to 2100.
- **Analysis Period**: Averaging statistics over 30-year intervals for a smooth mean time series spanning 1986 to 2085.
- **Resolution**: 0.1° in latitude and longitude, approximately equivalent to 11 kilometers at the equator.

### Exposure Data

Following the hazard identification, we will develop the exposure dataset. Exposure refers to assets, individuals, infrastructure, and other elements within a specified area that are susceptible to hazards, inclusive of their geographical coordinates, values, and additional relevant details.

**CLIMADA's Ready-to-Use Data: LitPop**
- **Description**: LitPop is a pre-existing exposure dataset provided by CLIMADA, integrating estimates of asset value, economic activity, or population based on nightlight intensity and population count data.
- **More Information**: [Click here for more information](https://www.climateadaptation.cc/tool/climada)

**External Exposure Data: Custom Dataset Creation**
- **Objective**: To create a custom exposure dataset, we will download and merge population data for Switzerland and Germany for the year 2020, offering data at 1 km resolution.
- **Steps to Create Exposure Data**:
  - **Data Source**: Access the Worldpop website (worldpop.org) or a specific data repository designated for population data.
  - **Data Selection**: Locate and download the population data for Switzerland and Germany at 1 km resolution for 2020.
  - **Format**: Ensure the data is in a format that aligns with your analytical or processing needs (e.g., raster file format).
 
## Environment Setup

To ensure consistency and reproducibility across different setups, we provide a Conda environment file, `climada_env_cds.yml`, that specifies all the dependencies needed to run the exercises and analyses presented in this project.

### Using the Provided Environment File

To recreate the project environment, follow these steps:

- **Ensure Conda is Installed**: Make sure you have Conda installed on your system. You can download Miniconda or Anaconda from their respective websites.

- **Download the Environment File**: The `climada_env_cds.yml` file is available in the root directory of this project. Download this file to your local machine.

- **Create the Environment**: Open your terminal (or Anaconda Prompt on Windows) and navigate to the directory where `climada_env_cds.yml` is located. Run the following command to create the Conda environment:
  ```bash
  conda env create -f climada_env_cds.yml

- **Activate the Environment**: Once the environment is successfully created, you can activate it using:
  ```bash
  conda activate climada_env_cds

- **Verify the Environment**: To ensure everything is set up correctly, you can list the installed packages using:
  ```bash
  conda list

This environment includes all necessary libraries and packages, such as CLIMADA, required data processing and analysis tools, and other dependencies specified within the project. By following these steps, you ensure that you have a consistent setup that mirrors the project's intended environment, facilitating a smoother workflow and analysis process.



