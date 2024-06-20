# <img src="src/TUIC.svg" alt="TUIC" width="40" height="40"> Sinkhole Cause Analysis
This project aimed to find important cause of road crisis and to utilize various open data to predict the happening of road crisis.

# Feature
- **Open data ETL** | Integrated 23 datasets from 4 opendata platform into different training features and analysis. The format of datasets includes csv, xml, shp, geojson and geopackage. A Selenium crawler to collect building cases information wich used to analyze the possibility that sinkhole case happened next to a building case.
- **Data preprocess** | Process spatial and time-series training features based on different geographic scale. (eg. Case location / Villages / 5m Hexes across Taipei City)
- **Model training** | Using XGBoost and LightGBM to find out important features.
- **Visualization** | Using Tableau, QGIS and also matplotlib to demonstrate important findings.

## Quick Start
1. install python 3.X and anaconda
2. create environment using anaconda and [Yaml file](https://github.com/Chu-c-git/Data-Science-Project_01_Roadcrisis_prevention_analysis/blob/main/practice_02_environment.yml).
   ```
   conda env create --file environment_name.yaml
   ```
3. Run related ipynb to make your own training data.

## Tools
| Tool | Description |
|---|---|
| Selenium | Automated web information crawling |
| Dask | Parallel processing for Python |
| Geopandas | Geographic data analysis in Python |
| XGBoost | Scalable and efficient gradient boosting machine learning |
| LightGBM | Gradient boosting machine learning with high performance |
| Optuna | Hyperparameter optimization library for machine learning |
| QGIS | Spatial visualization |