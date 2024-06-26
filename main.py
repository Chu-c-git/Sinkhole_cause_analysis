"""
main.py

This script orchestrates the entire workflow for downloading, preprocessing, and training a time 
series model using sinkhole data. It imports necessary libraries and modules, initializes the 
respective classes for downloading, preprocessing, and training, measures the time taken for each 
step, and logs the start and end times of each process along with the processing time.

Modules and Classes:
- os: Provides operating system-dependent functionality.
- xml.etree.ElementTree as ET: For parsing XML files.
- datetime from datetime: Used for date and time operations.
- geopandas as gpd: Handles geospatial data.
- function: Contains additional utility functions.
- download: Manages the download process of sinkhole data.
- preprocess_time_series: Handles preprocessing tasks such as data cleaning and transformation.
- train_time_series: Orchestrates the training of a time series model using XGBoost.

Workflow:
1. Downloading Data:
   - Initializes the download manager (`download.Download()`).
   - Measures and prints the start time of data downloading.
   - Runs the data downloading step (`download_manager.run()`).
   - Measures and prints the end time of data downloading.
   - Calculates and prints the processing time for data downloading.

2. Preprocessing Data:
   - Initializes the preprocessing manager (`preprocess_time_series.PreprocessTimeSeries()`).
   - Measures and prints the start time of data preprocessing.
   - Runs the data preprocessing step (`preprocess_time_series.run()`).
   - Measures and prints the end time of data preprocessing.
   - Calculates and prints the processing time for data preprocessing.

3. Model Training:
   - Initializes the model training manager (`train_time_series.TrainTimeSeries()`).
   - Measures and prints the start time of model training.
   - Runs the model training step (`train_time_series.run()`).
   - Measures and prints the end time of model training.
   - Calculates and prints the processing time for model training.
"""

import os
import xml.etree.ElementTree as ET
from datetime import datetime

import geopandas as gpd
from function import *
import download
import preprocess_time_series
import train_time_series

# Downloading
# Initialize with download
download_manager = download.Download()

# Measure and print the start time of data preprocessing
start_time = datetime.now()
print(f'Data Downloading Start at {start_time.strftime("%Y-%m-%d %H:%M:%S")}!\n')

# Run the data downloading step
download_manager.run()

# Measure and print the end time of data preprocessing
end_time = datetime.now()
print(f'\nData Downloading Finished at {end_time.strftime("%Y-%m-%d %H:%M:%S")}!')

# Calculate and print the processing time for data preprocessing
processing_time = end_time - start_time
print(f'\nProcessing Time on Data Downloading is {processing_time}!')


# Preprocessing
# Initialize with preprocess_time_series
preprocess_time_series = preprocess_time_series.PreprocessTimeSeries()

# Measure and print the start time of data preprocessing
start_time = datetime.now()
print(f'Data Preprocessing Start at {start_time.strftime("%Y-%m-%d %H:%M:%S")}!\n')

# Run the data downloading step
preprocess_time_series.run()

# Measure and print the end time of data preprocessing
end_time = datetime.now()
print(f'\nData Preprocessing Finished at {end_time.strftime("%Y-%m-%d %H:%M:%S")}!')

# Calculate and print the processing time for data preprocessing
processing_time = end_time - start_time
print(f'\nProcessing Time on Data Preprocessing is {processing_time}!')


# Model Training
# Initialize with train_time_series
train_time_series = train_time_series.TrainTimeSeries()

# Measure and print the start time of data preprocessing
start_time = datetime.now()
print(f'Model Training Started at {start_time.strftime("%Y-%m-%d %H:%M:%S")}!\n')

# Run the data downloading step
train_time_series.run()

# Measure and print the end time of data preprocessing
end_time = datetime.now()
print(f'\nModel Training Finished at {end_time.strftime("%Y-%m-%d %H:%M:%S")}!')

# Calculate and print the processing time for data preprocessing
processing_time = end_time - start_time
print(f'\nProcessing Time on Model Training is {processing_time}!')