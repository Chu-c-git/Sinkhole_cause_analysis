import json
import os
import subprocess
import time
import zipfile
from pathlib import Path

import fiona
import geopandas as gpd
import pandas as pd
import requests

def download_file(
    file_name,
    url,
    file_folder,
    is_proxy=False,
    is_verify=True,
    timeout: int = 60,
):
    """
    Download file from `url` to `{DATA_PATH}/{file_name}`.

    Args:
        file_name: str, file name
        url: str, file url
        is_proxy: bool, whether use proxy
        is_verify: bool, whether verify ssl
        timeout: int, request timeout
        file_folder: str, file folder path

    Returns: str, full file path

    Example:
    ``` python
    # Read GeoJSON
    # GeoJSON is a special format of JSON that represents geographical data
    # The extension of a geoJSON file can be .geojson or .json.
    import geopandas as gpd
    from utils.extract_stage import download_file_by_filename

    URL = "https://pwdgis.taipei/wg/opendata/I0201-5.geojson"
    FILE_FOLDER = "your_foler_path"
    FILE_NAME = "goose_sanctuary.geojson"
    FILE_ENCODING = "UTF-8"

    local_file = download_file_by_filename(FILE_NAME, URL, file_folder=FILE_FOLDER)
    gdata = gpd.read_file(local_file, encoding=FILE_ENCODING, driver="GeoJSON")
    print(gdata)
    ```
    ```
    output:
    Id     名稱            面積    類型  集水區  物理型  水文HY  濱水植  水質WQ  生物BI  MIWC2017                                           geometry
    0   3  雁鴨保護區  1.799444e+06  重要濕地  NaN  NaN   NaN  NaN   NaN   NaN       NaN  MULTIPOLYGON (((121.51075 25.02214, 121.51083 ...

    ```
    """
    full_file_path = f"{file_folder}/{file_name}"
    # download file
    try:
        with requests.get(
            url,
            stream=True,
            proxies=PROXIES if is_proxy else None,
            verify=is_verify,
            timeout=timeout,
        ) as r:
            r.raise_for_status()
            with open(full_file_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
        print(f"Downloaded {file_name} from {url}")
        return full_file_path
    except Exception as e:
        raise e
    
# fill missing level
def fill_missing_level(group):
    for i in range(len(group)):
        if pd.isna(group.iloc[i]['水位(m)']):
            # Get the water level values from the previous and next day
            prev_val = group.iloc[i-1]['水位(m)'] if i-1 >= 0 else None
            next_val = group.iloc[i+1]['水位(m)'] if i+1 < len(group) else None
            
            # Check if the previous and next day values are not missing
            if pd.notna(prev_val) and pd.notna(next_val):
                group.at[group.index[i], '水位(m)'] = (prev_val + next_val) / 2
            else:
                # Get the water level values from the previous two days and the next two days
                prev_vals = group.iloc[max(0, i-2):i]['水位(m)'].dropna().tolist()
                next_vals = group.iloc[i+1:i+3]['水位(m)'].dropna().tolist()
                combined_vals = prev_vals + next_vals
                if combined_vals:
                    group.at[group.index[i], '水位(m)'] = sum(combined_vals) / len(combined_vals)
    return group