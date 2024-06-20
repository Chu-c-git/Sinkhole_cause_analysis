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