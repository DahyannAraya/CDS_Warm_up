# Warm-Up Exercise - Integration of CLIMADA Impact Assessment with Copernicus Climate Change Service (C3S) data

**This initial exercise is designed as a framework to introduce users to Copernicus data and CLIMADA. Through this exercise, users will learn**:

- How to Search for Data: Guidelines on locating the required data sets.
- How to Access Data: Steps to gain access to data once found.
- Using Copernicus API with CLIMADA ETH: Instructions on utilizing the Copernicus API within the context of CLIMADA ETH for environmental data analysis.
- Data Exploration: Techniques for exploring the data to understand its structure and content.
- Introduction to CLIMADA ETH for Impact Analysis: An overview of how CLIMADA ETH can be used for assessing environmental impacts.
- Data Manipulation: How to process data to transform it into a hazard object suitable for analysis.
- Exploring the Hazard Object: Detailed guide on analyzing the hazard object to extract meaningful insights.
- Accessing Ready-to-Use Exposure Data (LitPop): Procedure to access LitPop, a ready-to-use exposure data set provided by CLIMADA.
- Accessing Custom Exposure Data: Instructions on how to search for and use your own exposure data for analysis.
- Developing an Impact Function: Steps to create a function to calculate the potential impact based on the data.
- Calculating the Impact: Methodology for computing the impact using the developed function.
- Exploring Impact Results: Techniques for analyzing and understanding the results of the impact assessment.
- Saving Analysis Results: Guidance on how to save the results of the analysis in both tabular and spatial formats for further use.

## Functions

Within the context of this project, specialized functions have been developed to assist users in exporting their findings, particularly the Annual Expected Impact (AEI). AEI refers to the anticipated impact in an average year for each exposure point, available in NetCDF or GeoTIFF formats. This feature is particularly beneficial for integrating CLIMADA's outcomes with other geospatial platforms or software for extended analysis. Moreover, these functionalities facilitate quicker and more straightforward visual communication of the results.

This exercise aims to equip users with the necessary skills to navigate, utilize, and analyze Copernicus data effectively within the CLIMADA ETH framework. By engaging with this guide, users will be prepared to conduct comprehensive impact analyses and communicate their findings efficiently.

## Data Description# Warm-Up Exercise - Integration of CLIMADA Impact Assessment with Copernicus Climate Change Service (C3S) data

**This initial exercise is designed as a framework to introduce users to Copernicus data and CLIMADA. Through this exercise, users will learn:**

- **How to Search for Data**: Guidelines on locating the required data sets.
- **How to Access Data**: Steps to gain access to data once found.
- **Using Copernicus API with CLIMADA ETH**: Instructions on utilizing the Copernicus API within the context of CLIMADA ETH for environmental data analysis.
- **Data Exploration**: Techniques for exploring the data to understand its structure and content.
- **Introduction to CLIMADA ETH for Impact Analysis**: An overview of how CLIMADA ETH can be used for assessing environmental impacts.
- **Data Manipulation**: How to process data to transform it into a hazard object suitable for analysis.
- **Exploring the Hazard Object**: Detailed guide on analyzing the hazard object to extract meaningful insights.
- **Accessing Ready-to-Use Exposure Data (LitPop)**: Procedure to access LitPop, a ready-to-use exposure data set provided by CLIMADA.
- **Accessing Custom Exposure Data**: Instructions on how to search for and use your own exposure data for analysis.
- **Developing an Impact Function**: Steps to create a function to calculate the potential impact based on the data.
- **Calculating the Impact**: Methodology for computing the impact using the developed function.
- **Exploring Impact Results**: Techniques for analyzing and understanding the results of the impact assessment.
- **Saving Analysis Results**: Guidance on how to save the results of the analysis in both tabular and spatial formats for further use.

## Functions

Within the context of this project, specialized functions have been developed to assist users in exporting their findings, particularly the Annual Expected Impact (AEI). AEI refers to the anticipated impact in an average year for each exposure point, available in NetCDF or GeoTIFF formats. This feature is particularly beneficial for integrating CLIMADA's outcomes with other geospatial platforms or software for extended analysis. Moreover, these functionalities facilitate quicker and more straightforward visual communication of the results.

This exercise aims to equip users with the necessary skills to navigate, utilize, and analyze Copernicus data effectively within the CLIMADA ETH framework. By engaging with this guide, users will be prepared to conduct comprehensive impact analyses and communicate their findings efficiently.

## Data Description

### Heatwave Hazard Analysis in Europe

**Objective**: In this first exercise, our goal is to analyze the heatwave hazard in Europe using the heatwave days dataset from Copernicus. This dataset is instrumental in understanding the potential heatwave impacts on the region, informed by climate change projections.

**Dataset Overview**:
- **Source**: Copernicus
- **Climate Scenarios**: RCP 8.5 and RCP 4.5, representing different greenhouse gas concentration trajectories.
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

Activate the Environment: Once the environment is successfully created, you can activate it using:
  ```bash
  conda activate climada_env_cds

Verify the Environment: To ensure everything is set up correctly, you can list the installed packages using:
  ```bash
  conda list

This environment includes all necessary libraries and packages, such as CLIMADA, required data processing and analysis tools, and other dependencies specified within the project. By following these steps, you ensure that you have a consistent setup that mirrors the project's intended environment, facilitating a smoother workflow and analysis process.**
