"""
download.py - Download data from the open data portal.

This script downloads data from the open data portal and saves it to the data/raw_data folder.
The download links are read from the open_data_download_url.csv file in the data folder. If the
csv file is not found, the script will print an error message and exit.

"""

import csv
import os
import time
import zipfile

import requests
from function import *
from tqdm import tqdm

# Set path
DOWNLOAD_FOLDER = 'data/raw_data'
CSV_FILE = 'data/open_data_download_url.csv'

# Function
## Create download folder
def create_download_folder():
    if not os.path.exists(DOWNLOAD_FOLDER):
        os.makedirs(DOWNLOAD_FOLDER)
    else:
        print('Folder exists.')

## Batch Download file
def download_file(file_name, url, download_folder, is_proxy=False, is_verify=True, timeout=60):
    file_path = os.path.join(download_folder, file_name)
    if os.path.exists(file_path):
        print(f'{file_name} already exists. Skipping download.')
        return

    try:
        with requests.get(url, stream=True, proxies=None if not is_proxy else {'http': 'http://your_proxy', 'https': 'https://your_proxy'}, verify=is_verify, timeout=timeout) as r:
            r.raise_for_status()
            with open(file_path, 'wb') as f:
                for chunk in tqdm(r.iter_content(chunk_size=8192), total=None, desc=f'Downloading {file_name}'):
                    f.write(chunk)
    except requests.exceptions.RequestException as e:
        print(f"Error downloading {file_name} from {url}: {e}")

## Check missing files
def check_missing_files(download_folder, expected_files):
    downloaded_files = set(os.listdir(download_folder))
    missing_files = [file for file in expected_files if file not in downloaded_files]
    if missing_files:
        print("Missing files:")
        for file in missing_files:
            print(file)
    else:
        print("All files have been downloaded.")

## Unzip files
def unzip_files(download_folder):
    for item in os.listdir(download_folder):
        if item.endswith('.zip'):
            file_path = os.path.join(download_folder, item)
            folder_name = os.path.splitext(item)[0]
            unzip_dir = os.path.join(download_folder, folder_name)
            os.makedirs(unzip_dir, exist_ok=True)
            with zipfile.ZipFile(file_path, 'r') as zip_ref:
                zip_ref.extractall(unzip_dir)
            print(f"{item} unzipped to {folder_name}/")

## Main
def main():
    create_download_folder()

    try:
        url_file = pd.read_csv(CSV_FILE, encoding='utf-8')
    except FileNotFoundError:
        print("CSV file not found.")
        return

    auto_df = url_file[url_file['download_type'] == 'Auto']
    manual_df = url_file[url_file['download_type'] == 'Manually']
    request_df = url_file[url_file['download_type'] == 'Request']

    for idx, row in auto_df.iterrows():
        download_file(row['file_name'], row['download_link'], DOWNLOAD_FOLDER)

    unzip_files(DOWNLOAD_FOLDER)

    print("== Auto_download_file ==")
    check_missing_files(DOWNLOAD_FOLDER, auto_df['file_name'])

    print("== Manually_download_file ==")
    check_missing_files(DOWNLOAD_FOLDER, manual_df['file_name'])

    print("== Requested_download_file ==")
    check_missing_files(DOWNLOAD_FOLDER, request_df['file_name'])

if __name__ == "__main__":
    main()