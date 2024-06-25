"""
download.py - Download data from the open data portal.

This script downloads data from the open data portal and saves it to the data/raw_data folder.
The download links are read from the open_data_download_url.csv file in the data folder. If the
csv file is not found, the script will print an error message and exit.

"""

import csv
import os
import zipfile
from datetime import datetime

import pandas as pd
import requests
from tqdm import tqdm


class Download:
    def __init__(
        self,
        download_folder="data/raw_data",
        csv_file="data/open_data_download_url.csv",
    ):
        self.download_folder = download_folder
        self.csv_file = csv_file
        self.create_download_folder()

    def create_download_folder(self):
        if not os.path.exists(self.download_folder):
            os.makedirs(self.download_folder)
        else:
            print("Folder exists.")

    def download_file(self, file_name, url, is_proxy=False, is_verify=True, timeout=60):
        file_path = os.path.join(self.download_folder, file_name)
        if os.path.exists(file_path):
            print(f"{file_name} already exists. Skipping download.")
            return

        try:
            with requests.get(
                url,
                stream=True,
                proxies=(
                    None
                    if not is_proxy
                    else {"http": "http://your_proxy", "https": "https://your_proxy"}
                ),
                verify=is_verify,
                timeout=timeout,
            ) as r:
                r.raise_for_status()
                with open(file_path, "wb") as f:
                    for chunk in tqdm(
                        r.iter_content(chunk_size=8192),
                        total=None,
                        desc=f"Downloading {file_name}",
                    ):
                        f.write(chunk)
        except requests.exceptions.RequestException as e:
            print(f"Error downloading {file_name} from {url}: {e}")

    def check_missing_files(self, expected_files):
        downloaded_files = set(os.listdir(self.download_folder))
        missing_files = [
            file for file in expected_files if file not in downloaded_files
        ]
        if missing_files:
            print("Missing files:")
            for file in missing_files:
                print(file)
        else:
            print("All files have been downloaded.")

    def unzip_files(self):
        for item in os.listdir(self.download_folder):
            if item.endswith(".zip"):
                file_path = os.path.join(self.download_folder, item)
                folder_name = os.path.splitext(item)[0]
                unzip_dir = os.path.join(self.download_folder, folder_name)
                os.makedirs(unzip_dir, exist_ok=True)
                with zipfile.ZipFile(file_path, "r") as zip_ref:
                    zip_ref.extractall(unzip_dir)
                print(f"{item} unzipped to {folder_name}/")

    def run(self):
        start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"Data Downloading Start at {start_time}!\n")

        try:
            url_file = pd.read_csv(self.csv_file, encoding="utf-8")
        except FileNotFoundError:
            print("CSV file not found.")
            return

        auto_df = url_file[url_file["download_type"] == "Auto"]
        manual_df = url_file[url_file["download_type"] == "Manually"]
        request_df = url_file[url_file["download_type"] == "Request"]

        for idx, row in auto_df.iterrows():
            self.download_file(row["file_name"], row["download_link"])

        self.unzip_files()

        print("== Auto_download_file ==")
        self.check_missing_files(auto_df["file_name"])

        print("== Manually_download_file ==")
        self.check_missing_files(manual_df["file_name"])

        print("== Requested_download_file ==")
        self.check_missing_files(request_df["file_name"])

        end_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"\nData Downloading Complete at {end_time}!")


# 如果需要直接執行這個腳本，可以加上以下程式碼
if __name__ == "__main__":
    download_manager = Download()
    download_manager.run()
