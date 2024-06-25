# <img src="src/TUIC.svg" alt="TUIC" width="30" height="30"> Sinkhole Cause Analysis
This project aimed to find important cause of road crisis and to utilize various open data to predict the happening of road crisis.

# Feature
- **Open data ETL** | Integrated 23 datasets from 4 opendata platform into different training features and analysis. The format of datasets includes csv, xml, shp, geojson and geopackage. 
- **Data preprocess** | Process spatial and time-series training features based on different geographic scale. (eg. Case location / 5m Hexes across Taipei City)
- **Model training** | Using XGBoost and LightGBM to find out important features.
- **Visualization** | Using Tableau, QGIS and also matplotlib to demonstrate important findings.

# License
This project is licensed under the Apache 2.0 License. For more details, please refer to the LICENSE file.

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
| Dask | Parallel processing for Python |
| Geopandas | Geographic data analysis in Python |
| XGBoost | Scalable and efficient gradient boosting machine learning |
| LightGBM | Gradient boosting machine learning with high performance |
| Optuna | Hyperparameter optimization library for machine learning |
| QGIS | Spatial visualization |