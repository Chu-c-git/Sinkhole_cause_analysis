# <img src="src/TUIC.svg" alt="TUIC" width="30" height="30"> Sinkhole Cause Analysis
## Description
This project is a collaboration between the New Construction Office, Public Works Department, Taipei City Government and the Taipei Urban Intelligence Center. By utilizing open data, the project aims to analyze and discuss potential causes of hazardous sinkholes, allowing for early prevention.

Using machine learning models, the project conducts time series analysis on a monthly basis to identify potential causes of sinkholes. Various open data sources, such as tidal information, groundwater levels, and earthquake occurrences, are integrated into the research.

The final results from the machine learning models can predict the number of sinkholes. The model's findings on the importance of different causes will provide feedback for practical discussions and strategies.

## Dataset
### Road Repair Case Data
Records of road repair-related cases, including details reported by the public, case coordinates, and the status of construction and dispatch.

### Water Pipeline Repair Case Data
Records of water pipeline (such as tap water pipelines) repair cases, including the construction sections, report dates, and repair purposes.

### Open data
For the usage of open data, you can refer to the list [here.](https://github.com/Chu-c-git/Sinkhole_cause_analysis/blob/main/data/open_data_download_url.csv) The data includes both time-series and spatial analysis. However, for now, the repository only releases the time-series module.

## Quick Start
1. install python 3.X and anaconda
2. create environment using anaconda and [Yaml file](https://github.com/Chu-c-git/Sinkhole_cause_analysis/blob/main/environment.yaml).
   Run the following command to setup the conda enviroment.
   ``` bash
   conda env create -f /path/to/environment.yml
   ```
3. Run related ipynb to make your own training data.

## Feature
- **Open data ETL** | Integrated 23 datasets from 4 opendata platform into different training features and analysis. The format of datasets includes csv, xml, shp, geojson and geopackage. 
- **Data preprocess** | Process spatial and time-series training features based on different geographic scale. (eg. Case location / 5m Hexes across Taipei City)
- **Model training** | Using XGBoost to find out important features.
- **Visualization** | Using Tableau, QGIS and also matplotlib to demonstrate important findings.

## License
This project is licensed under the Apache 2.0 License. For more details, please refer to the LICENSE file.

## Acknowledgements
The authors thank [**New Construction Office, Public Works Department, Taipei City Government**](https://nco.gov.taipei/Default.aspx) for their invaluable collaboration throughout the project, from its initiation, through ongoing communication, to its successful execution. They played a crucial role in this endeavor. We would also like to express our gratitude to [**Road Excavation Administration Center, Public Works Department, Taipei City Government**](https://dig.taipei/Tpdig/), and [**Water Resources Agency, MOEA**](https://gweb.wra.gov.tw/Hydroinfo/) for providing the data that enabled us to complete this project.

# References
1. [臺北市資料大平台](https://data.taipei/)
2. [政府資料開放平台](https://data.gov.tw/)
3. [經濟部水利署水文資訊網](https://gweb.wra.gov.tw/Hydroinfo/)
4. [XGBoost](https://xgboost.readthedocs.io/en/stable/)
5. [Optuna](https://optuna.org/)