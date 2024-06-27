import json
import os
import subprocess
import time
import xml.etree.ElementTree as ET
import zipfile
from datetime import datetime, timedelta
from pathlib import Path

import dask.dataframe as dd
import dask_geopandas as dgpd
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import requests
import seaborn as sns
import tqdm
from dask.diagnostics import ProgressBar
from matplotlib.patches import Polygon
from pyproj import Proj
from shapely import LineString
from shapely.geometry import LineString, Point, Polygon
from shapely.wkt import loads
from tqdm import tqdm


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
        if pd.isna(group.iloc[i]["水位(m)"]):
            # Get the water level values from the previous and next day
            prev_val = group.iloc[i - 1]["水位(m)"] if i - 1 >= 0 else None
            next_val = group.iloc[i + 1]["水位(m)"] if i + 1 < len(group) else None

            # Check if the previous and next day values are not missing
            if pd.notna(prev_val) and pd.notna(next_val):
                group.at[group.index[i], "水位(m)"] = (prev_val + next_val) / 2
            else:
                # Get the water level values from the previous two days and the next two days
                prev_vals = group.iloc[max(0, i - 2) : i]["水位(m)"].dropna().tolist()
                next_vals = group.iloc[i + 1 : i + 3]["水位(m)"].dropna().tolist()
                combined_vals = prev_vals + next_vals
                if combined_vals:
                    group.at[group.index[i], "水位(m)"] = sum(combined_vals) / len(
                        combined_vals
                    )
    return group


def process_soil_liquid(hex_5m_sample, soil_liquid, epsg_code=3826):
    """
    Process soil_liquid GeoDataFrame:
    1. Convert CRS to EPSG:3826.
    2. Perform spatial join with hex_5m_sample based on 'within'.
    3. Drop unnecessary columns.
    4. Rename the 'class' column to 'soil_liquid_class'.
    5. Replace the soil_liquid_class with 0, 1, 2.
    6. Drop duplicates.

    Parameters:
    - hex_5m_sample (GeoDataFrame): GeoDataFrame containing hexagons.
    - soil_liquid (GeoDataFrame): GeoDataFrame containing soil_liquid data.
    - epsg_code (int): EPSG code for the target CRS (default is 3826).

    Returns:
    - GeoDataFrame: Resulting GeoDataFrame after processing.
    """
    # Convert CRS to EPSG:3826
    soil_liquid = soil_liquid.to_crs(epsg=epsg_code)

    # Spatial join based on 'within'
    hex_5m_rd_sample_s = gpd.sjoin(
        hex_5m_sample, soil_liquid, how="inner", predicate="intersects"
    )

    # Drop unnecessary columns
    drop_col = ["index_right"]
    hex_5m_rd_sample_s.drop(columns=drop_col, inplace=True)

    # Rename the 'class' column
    rename_col = {"class": "soil_liquid_class"}
    hex_5m_rd_sample_s.rename(columns=rename_col, inplace=True)

    # Replace the soil_liquid_class
    soil_replace = {"1": 2, "2": 1, "3": 0}
    hex_5m_rd_sample_s["soil_liquid_class"] = hex_5m_rd_sample_s[
        "soil_liquid_class"
    ].replace(soil_replace)

    # Drop duplicates
    hex_5m_rd_sample_s.drop_duplicates(subset=["id"], inplace=True)

    return hex_5m_rd_sample_s


def road_properties_process(above_8m_road, under_8m_road, hex_5m_rd_sample_s):
    """
    Process road properties:
    1. Drop useless columns & rename columns for above_8m_road and under_8m_road.
    2. Concatenate above_8m_road and under_8m_road.
    3. Perform spatial join with hex_5m_rd_sample_s based on 'intersects'.
    4. Drop unnecessary columns.
    5. Drop duplicates.

    Parameters:
    - above_8m_road (GeoDataFrame): GeoDataFrame containing road properties for roads above 8m.
    - under_8m_road (GeoDataFrame): GeoDataFrame containing road properties for roads under 8m.
    - hex_5m_rd_sample_s (GeoDataFrame): GeoDataFrame containing hexagons within the Taipei roads with soil_liquid.

    Returns:
    - GeoDataFrame: Resulting GeoDataFrame after processing.
    """
    # Process above_8m_road
    keep_col_a = ["RoadWidth", "路名", "Road_ID", "geometry"]
    rename_col = {"RoadWidth": "width", "路名": "road_name", "Road_ID": "road_id"}
    above_8m_road = above_8m_road[keep_col_a].copy()
    above_8m_road.rename(columns=rename_col, inplace=True)

    # Process under_8m_road
    keep_col_u = ["ROADID", "ROADNAME", "WIDTH", "geometry"]
    rename_col = {"WIDTH": "width", "ROADNAME": "road_name", "ROADID": "road_id"}
    under_8m_road = under_8m_road[keep_col_u].copy()
    under_8m_road.rename(columns=rename_col, inplace=True)

    # Concatenate above_8m_road & under_8m_road
    road_all = pd.concat([above_8m_road, under_8m_road], ignore_index=True)

    # Perform spatial join based on 'intersects'
    hex_5m_rd_sample_r = gpd.sjoin(
        hex_5m_rd_sample_s, road_all, how="inner", predicate="intersects"
    )

    # Drop unnecessary columns
    drop_col = ["index_right"]
    hex_5m_rd_sample_r.drop(columns=drop_col, inplace=True)

    # Drop duplicates
    hex_5m_rd_sample_r.drop_duplicates(subset=["id"], inplace=True)

    return hex_5m_rd_sample_r


def calculate_pipeline_count(hexagon_df, pipeline_gdf, count_column="sp_count"):
    """
    計算每個 hexagon 中管線的數量並更新 DataFrame。

    Parameters:
    - hexagon_df (geopandas.GeoDataFrame): 包含 hexagon 的 GeoDataFrame。
    - pipeline_gdf (geopandas.GeoDataFrame): 包含管線的 GeoDataFrame。
    - count_column (str): 保存管線數量的欄位名稱。

    Returns:
    - geopandas.GeoDataFrame: 更新後的 GeoDataFrame。
    """
    hexagon_df_copy = hexagon_df.copy()

    for idx, row in tqdm(
        hexagon_df_copy.iterrows(), total=len(hexagon_df_copy), desc=count_column
    ):
        buffer_zone = row["geometry"]
        hexagon_df_copy.at[idx, count_column] = pipeline_gdf.intersects(
            buffer_zone
        ).sum()

    return hexagon_df_copy


def summarize_pipe_counts(hex_5m_rd_sample_cn):
    """
    Summarize pipeline counts for water and water-related systems.

    Parameters:
    - hex_5m_rd_sample_cn (DataFrame): DataFrame containing pipeline counts for each hexagon.

    Returns:
    - DataFrame: Resulting DataFrame after summarizing pipeline counts.
    """
    # 自來水管線數量彙整
    hex_5m_rd_sample_cn["wp_count"] = hex_5m_rd_sample_cn[
        ["wp_01_count", "wp_02_count", "wp_03_count", "wp_04_count"]
    ].sum(axis=1)
    drop_col_wp = ["wp_01_count", "wp_02_count", "wp_03_count", "wp_04_count"]
    hex_5m_rd_sample_cn.drop(columns=drop_col_wp, inplace=True)

    # 水系管線數量彙整
    hex_5m_rd_sample_cn["pipe_count"] = hex_5m_rd_sample_cn[
        ["sp_count", "rp_count", "rd_count", "wp_count", "cn_count"]
    ].sum(axis=1)

    return hex_5m_rd_sample_cn


def create_buffered_gdf(
    df, lon_col, lat_col, buffer_distance=5, geometry_col="geometry"
):
    """
    創建包含緩衝區的 GeoDataFrame。

    Parameters:
    - df (pandas.DataFrame): 包含案件經緯度的 DataFrame。
    - lon_col (str): 經度欄位名稱。
    - lat_col (str): 緯度欄位名稱。
    - buffer_distance (float): 緩衝區的距離。
    - geometry_col (str): 存儲緩衝區的欄位名稱。

    Returns:
    - geopandas.GeoDataFrame: 包含緩衝區的 GeoDataFrame。
    """
    df = df.copy()
    # 將經緯度轉換為 Point
    df.loc[:, "centroid"] = df.apply(lambda x: Point(x[lon_col], x[lat_col]), axis=1)

    # 創建緩衝區
    df.loc[:, geometry_col] = df.apply(
        lambda x: x.centroid.buffer(buffer_distance), axis=1
    )

    # 將 DataFrame 轉換為 GeoDataFrame
    geo_df = gpd.GeoDataFrame(df, geometry=geometry_col)
    geo_df = geo_df.drop(columns=["centroid"])

    return geo_df


def calculate_case_count(hexagon_df, case_gdf, count_column="sp_count"):
    """
    計算每個 hexagon 中管線的數量並更新 DataFrame。

    Parameters:
    - hexagon_df (geopandas.GeoDataFrame): 包含 hexagon 的 GeoDataFrame。
    - case_gdf (geopandas.GeoDataFrame): 包含案件的 GeoDataFrame。
    - count_column (str): 保存管線數量的欄位名稱。

    Returns:
    - geopandas.GeoDataFrame: 更新後的 GeoDataFrame。
    """
    hexagon_df_copy = hexagon_df.copy()
    # Setting crs and columns to drop
    hexagon_df_copy.crs = "EPSG:3826"
    case_gdf.crs = "EPSG:3826"
    keep_col = ["geometry"]
    drop_col = [col for col in case_gdf.columns.tolist() if col not in keep_col]

    # Spatial join
    hexagon_df_copy = gpd.sjoin(
        hexagon_df_copy, case_gdf, how="left", predicate="intersects"
    )
    hexagon_df_copy = hexagon_df_copy.drop(columns=drop_col)

    # Counting
    hexagon_df_copy["index_right"] = (
        hexagon_df_copy["index_right"].fillna(0).astype(int)
    )
    count_result = hexagon_df_copy.groupby("id").size()
    hexagon_df_copy[count_column] = hexagon_df_copy["id"].map(count_result)

    # Fill 0 intersection
    hexagon_df_copy[count_column] = hexagon_df_copy[count_column].where(
        hexagon_df_copy["index_right"] != 0, 0
    )
    right_col = [col for col in hexagon_df_copy.columns if col.endswith("right")]
    hexagon_df_copy = hexagon_df_copy.drop(columns=right_col)
    hexagon_df_copy = hexagon_df_copy.rename(columns={"geometry_left": "geometry"})
    hexagon_df_copy = hexagon_df_copy.drop_duplicates(subset=["id"])
    hexagon_df_copy = gpd.GeoDataFrame(hexagon_df_copy, geometry="geometry")

    # print(f"{count_column} is done!")
    return hexagon_df_copy


def calculate_case_count_v2(
    hexagon_df, case_gdf, count_column="sp_count", groupby="id"
):
    """
    計算每個 hexagon 中管線的數量並更新 DataFrame。

    Parameters:
    - hexagon_df (geopandas.GeoDataFrame): 包含 hexagon 的 GeoDataFrame。
    - case_gdf (geopandas.GeoDataFrame): 包含案件的 GeoDataFrame。
    - count_column (str): 保存管線數量的欄位名稱。

    Returns:
    - geopandas.GeoDataFrame: 更新後的 GeoDataFrame。
    """
    hexagon_df_copy = hexagon_df.copy()
    # Setting crs and columns to drop
    hexagon_df_copy.crs = "EPSG:3826"
    case_gdf.crs = "EPSG:3826"
    keep_col = ["geometry"]
    drop_col = [col for col in case_gdf.columns.tolist() if col not in keep_col]

    # Spatial join
    hexagon_df_copy = gpd.sjoin(
        hexagon_df_copy, case_gdf, how="left", predicate="intersects"
    )
    hexagon_df_copy = hexagon_df_copy.drop(columns=drop_col)

    # Counting
    hexagon_df_copy["index_right"] = (
        hexagon_df_copy["index_right"].fillna(0).astype(int)
    )
    count_result = hexagon_df_copy.groupby(groupby).size()
    hexagon_df_copy[count_column] = hexagon_df_copy[groupby].map(count_result)

    # Fill 0 intersection
    hexagon_df_copy[count_column] = hexagon_df_copy[count_column].where(
        hexagon_df_copy["index_right"] != 0, 0
    )
    right_col = [col for col in hexagon_df_copy.columns if col.endswith("right")]
    hexagon_df_copy = hexagon_df_copy.drop(columns=right_col)
    hexagon_df_copy = hexagon_df_copy.rename(columns={"geometry_left": "geometry"})
    hexagon_df_copy = hexagon_df_copy.drop_duplicates(subset=[groupby])
    hexagon_df_copy = gpd.GeoDataFrame(hexagon_df_copy, geometry="geometry")

    # print(f"{count_column} is done!")
    return hexagon_df_copy


def calculate_case_during_period(
    hexagon_df, case_gdf, col_name, date, time_window=7, count_column="sp_count"
):
    """
    計算每個 hexagon 中符合時間區間的案件數量並更新 DataFrame。

    Parameters:
    - hexagon_df (geopandas.GeoDataFrame): 包含 hexagon 的 GeoDataFrame。
    - case_gdf (geopandas.GeoDataFrame): 包含案件的 GeoDataFrame。
    - date (datetime.datetime): 預測日的日期。
    - time_window (int): 時間區間的長度。
    - count_column (str): 保存管線數量的欄位名稱。

    Returns:
    - geopandas.GeoDataFrame: 更新後的 GeoDataFrame。
    """
    window_start = date - timedelta(days=time_window + 1)
    window_end = date - timedelta(days=1)
    # case_gdf['查報日期'] = case_gdf['查報日期'].dt.date
    mask = (case_gdf[col_name] >= window_start) & (case_gdf[col_name] <= window_end)
    incidents_in_window = case_gdf[mask].copy()
    incidents_in_window.set_crs(epsg=3826, inplace=True)

    # Setting crs and columns to drop
    hexagon_df = hexagon_df.set_crs(epsg=3826)
    keep_col = ["geometry"]
    drop_col = [
        col for col in incidents_in_window.columns.tolist() if col not in keep_col
    ]

    # Spatial join
    hexagon_df_copy = gpd.sjoin(
        hexagon_df, incidents_in_window, how="left", predicate="intersects"
    )
    hexagon_df_copy = hexagon_df_copy.drop(columns=drop_col)

    # Counting
    hexagon_df_copy["index_right"] = (
        hexagon_df_copy["index_right"].fillna(0).astype(int)
    )
    count_result = hexagon_df_copy.groupby("id").size()
    hexagon_df_copy[count_column] = hexagon_df_copy["id"].map(count_result)

    # Fill 0 intersection
    hexagon_df_copy.loc[:, count_column] = hexagon_df_copy[count_column].where(
        hexagon_df_copy["index_right"] != 0, 0
    )
    right_col = [col for col in hexagon_df_copy.columns if col.endswith("right")]
    hexagon_df_copy = hexagon_df_copy.drop(columns=right_col)
    hexagon_df_copy = hexagon_df_copy.rename(columns={"geometry_left": "geometry"})
    hexagon_df_copy = hexagon_df_copy.drop_duplicates(subset=["id"])
    hexagon_df_copy = gpd.GeoDataFrame(hexagon_df_copy, geometry="geometry")

    # print(f"{count_column} is done!")
    return hexagon_df_copy


def calculate_case_during_period_boolean(
    hexagon_df, case_gdf, col_name, date, time_window=7, count_column="sp_count"
):
    """
    計算每個 hexagon 中符合時間區間的案件數量並更新 DataFrame。

    Parameters:
    - hexagon_df (geopandas.GeoDataFrame): 包含 hexagon 的 GeoDataFrame。
    - case_gdf (geopandas.GeoDataFrame): 包含案件的 GeoDataFrame。
    - date (datetime.datetime): 預測日的日期。
    - time_window (int): 時間區間的長度。
    - count_column (str): 保存管線數量的欄位名稱。

    Returns:
    - geopandas.GeoDataFrame: 更新後的 GeoDataFrame。
    """
    window_start = date - timedelta(days=time_window + 1)
    window_end = date - timedelta(days=1)
    # case_gdf['查報日期'] = case_gdf['查報日期'].dt.date
    mask = (case_gdf[col_name] >= window_start) & (case_gdf[col_name] <= window_end)
    incidents_in_window = case_gdf[mask].copy()
    incidents_in_window.set_crs(epsg=3826, inplace=True)

    # Setting crs and columns to drop
    hexagon_df = hexagon_df.set_crs(epsg=3826)
    keep_col = ["geometry"]
    drop_col = [
        col for col in incidents_in_window.columns.tolist() if col not in keep_col
    ]

    # Spatial join
    hexagon_df_copy = gpd.sjoin(
        hexagon_df, incidents_in_window, how="left", predicate="intersects"
    )
    hexagon_df_copy = hexagon_df_copy.drop(columns=drop_col)

    # Counting
    hexagon_df_copy["index_right"] = (
        hexagon_df_copy["index_right"].fillna(0).astype(int)
    )
    count_result = hexagon_df_copy.groupby("id").size()
    hexagon_df_copy[count_column] = hexagon_df_copy["id"].map(count_result)

    # Fill 0 intersection (boolean)
    hexagon_df_copy.loc[:, count_column] = hexagon_df_copy[count_column].where(
        hexagon_df_copy["index_right"] != 0, 0
    )
    right_col = [col for col in hexagon_df_copy.columns if col.endswith("right")]
    hexagon_df_copy = hexagon_df_copy.drop(columns=right_col)
    hexagon_df_copy = hexagon_df_copy.rename(columns={"geometry_left": "geometry"})
    hexagon_df_copy = hexagon_df_copy.drop_duplicates(subset=["id"])
    hexagon_df_copy[count_column] = hexagon_df_copy[count_column].apply(
        lambda x: 1 if x > 0 else 0
    )
    hexagon_df_copy = gpd.GeoDataFrame(hexagon_df_copy, geometry="geometry")

    # print(f"{count_column} is done!")
    return hexagon_df_copy


def calculate_case_on_date(
    hexagon_df, case_gdf, col_name, date, count_column="sp_count"
):
    """
    計算每個 hexagon 中符合時間區間的案件數量並更新 DataFrame。

    Parameters:
    - hexagon_df (geopandas.GeoDataFrame): 包含 hexagon 的 GeoDataFrame。
    - case_gdf (geopandas.GeoDataFrame): 包含案件的 GeoDataFrame。
    - date (datetime.datetime): 預測日的日期。
    - count_column (str): 保存管線數量的欄位名稱。

    Returns:
    - geopandas.GeoDataFrame: 更新後的 GeoDataFrame。
    """
    # case_gdf['查報日期'] = case_gdf['查報日期'].dt.date
    date = pd.to_datetime(date).date()
    mask = case_gdf[col_name] == date
    incidents_in_window = case_gdf[mask].copy()
    incidents_in_window.set_crs(epsg=3826, inplace=True)

    # Setting crs and columns to drop
    hexagon_df = hexagon_df.set_crs(epsg=3826)
    keep_col = ["geometry"]
    drop_col = [
        col for col in incidents_in_window.columns.tolist() if col not in keep_col
    ]

    # Spatial join
    hexagon_df_copy = gpd.sjoin(
        hexagon_df, incidents_in_window, how="left", predicate="intersects"
    )
    hexagon_df_copy = hexagon_df_copy.drop(columns=drop_col)

    # Counting
    hexagon_df_copy["index_right"] = (
        hexagon_df_copy["index_right"].fillna(0).astype(int)
    )
    count_result = hexagon_df_copy.groupby("id").size()
    hexagon_df_copy[count_column] = hexagon_df_copy["id"].map(count_result)

    # Fill 0 intersection
    hexagon_df_copy.loc[:, count_column] = hexagon_df_copy[count_column].where(
        hexagon_df_copy["index_right"] != 0, 0
    )
    right_col = [col for col in hexagon_df_copy.columns if col.endswith("right")]
    hexagon_df_copy = hexagon_df_copy.drop(columns=right_col)
    hexagon_df_copy = hexagon_df_copy.rename(columns={"geometry_left": "geometry"})
    hexagon_df_copy = hexagon_df_copy.drop_duplicates(subset=["id"])
    hexagon_df_copy = gpd.GeoDataFrame(hexagon_df_copy, geometry="geometry")

    # print(f"{count_column} is done!")
    return hexagon_df_copy


def rocdate_transfer_to_time(column):
    """
    将民國年日期轉換為時間序列

    Parameters:
    - column: pandas.Series，包含民國年日期的列

    Returns:
    - pandas.Series，轉換後的時間序列
    """
    # 將列轉換為字符串
    column = column.astype(str)

    # 拆分民國年與時間，並添加條件檢查，確保 roc_year 不是空字符串
    roc_year = column.str[:3]
    mask = roc_year != ""
    roc_year = roc_year[mask]

    if not roc_year.empty:
        # 民國轉西元
        month = column.str[3:5]
        day = column.str[5:]
        ad_year = roc_year.astype(int) + 1911
        ad_date = ad_year.astype(str).str.cat([month, day], sep="/")
        # 轉換為時間格式
        new_date = pd.to_datetime(ad_date, format="%Y/%m/%d").dt.date
        column = new_date

    return column


def find_average_rainfall(
    time_data, time_start, time_end, column="daily_precipitation"
):
    # Create mask based on time range
    mask = (time_data["date"] >= time_start) & (time_data["date"] < time_end)

    # Apply mask to create a subset of data
    time_data_mask = time_data.loc[mask]

    # Calculate the average rainfall
    average_rainfall = time_data_mask[column].mean().round(3)

    return average_rainfall


def process_rainfall_data(
    time_data, time_start, time_pred, hex_5m_rd_sample_cn, column="precipitation"
):
    # Calculate average rainfall using the provided function
    average_rainfall = find_average_rainfall(
        time_data, time_start, time_pred, column="daily_precipitation"
    )

    # Update the "precipitation" column in hex_5m_rd_sample_cn
    hex_5m_rd_sample_cn[column] = average_rainfall

    return hex_5m_rd_sample_cn


def process_sum_data(
    time_data, time_start, time_pred, hex_5m_rd_sample_cn, column="precipitation"
):
    # Create mask based on time range
    mask = (time_data["date"] >= time_start) & (time_data["date"] < time_pred)

    # Apply mask to create a subset of data
    time_data_mask = time_data.loc[mask]

    # Calculate the average rainfall
    sum_case_count = time_data_mask[column].sum()

    # Update the "precipitation" column in hex_5m_rd_sample_cn
    hex_5m_rd_sample_cn[column] = sum_case_count

    return hex_5m_rd_sample_cn


def process_mean_data(
    time_data, time_start, time_pred, hex_5m_rd_sample_cn, column="precipitation"
):
    # Create mask based on time range
    mask = (time_data["date"] >= time_start) & (time_data["date"] < time_pred)

    # Apply mask to create a subset of data
    time_data_mask = time_data.loc[mask]
    time_data_mask.loc[:, column] = pd.to_numeric(time_data_mask[column]).copy()

    # Calculate the average rainfall
    sum_case_count = time_data_mask[column].mean().round(3)

    # Update the "precipitation" column in hex_5m_rd_sample_cn
    hex_5m_rd_sample_cn[column] = sum_case_count

    return hex_5m_rd_sample_cn


def extract_lng(centroid):
    if isinstance(centroid, Point):
        return centroid.x
    else:
        return loads(centroid).x


def extract_lat(centroid):
    if isinstance(centroid, Point):
        return centroid.y
    else:
        return loads(centroid).y


# 列印沒有namespace的樹狀結構
def print_xml_tree_without_namespace(element, indent=0):
    # 取得不包含命名空間的 tag 名稱
    tag_name = element.tag.split("}")[1] if "}" in element.tag else element.tag

    # 顯示當前節點
    print("  " * indent + f"{tag_name}: {element.text}")

    # 遞迴顯示子節點
    for child in element:
        print_xml_tree_without_namespace(child, indent + 1)


# 列印有namespace的樹狀結構
def print_xml_tree_with_namespace(element, indent=0):
    # 顯示當前節點
    print("  " * indent + f"{element.tag}: {element.text}")

    # 遞迴顯示子節點
    for child in element:
        print_xml_tree_with_namespace(child, indent + 1)


def data_to_linestring(pos_list_text):

    # 將 posList 字串分割為座標列表
    coordinates = [float(coord) for coord in pos_list_text.split()]

    # 將座標列表轉換成符合需求的格式
    input_str = [
        (coordinates[i], coordinates[i + 1], coordinates[i + 2])
        for i in range(0, len(coordinates), 3)
    ]

    # 使用 eval 将字符串转换为 Python 对象
    coord = eval(str(input_str))

    # 確保 coords 中包含至少兩個點
    if len(coord) < 2:
        raise ValueError("Coordinate array must contain at least two points.")

    # 构建 LineString 对象
    line = LineString(coord)
    return str(line)


def remove_z_from_linestring(wkt_with_z):
    # 使用 Shapely 的 loads 方法解析 WKT 格式的几何对象
    line_with_z = loads(wkt_with_z)
    line_with_z = line_with_z.astype(str)

    # 通过访问坐标的 xy 属性创建新的 LineString，去除 Z 值
    coords_without_z = list(zip(line_with_z.xy[0], line_with_z.xy[1]))
    line_without_z = LineString(coords_without_z)

    return line_without_z


def plot_confusion_matrix(actual_val, pred_val, title=None):
    confusion_matrix = pd.crosstab(
        actual_val, pred_val, rownames=["Actual"], colnames=["Predicted"]
    )

    plot = sns.heatmap(confusion_matrix, annot=True, fmt=",.0f")

    if title is None:
        pass
    else:
        plot.set_title(title)

    # plt.show()


def plot_roc_curve(model, X_test, y_test):
    """
    繪製模型的 ROC 曲線並計算 ROC AUC 分數。

    Parameters:
    - model: 已經訓練好的分類模型，例如 XGBoost。
    - X_test: 測試集的特徵。
    - y_test: 測試集的目標變數。

    Returns:
    - None
    """
    # 使用模型預測測試集概率
    y_proba = model.predict_proba(X_test)[:, 1]

    # 計算 ROC 曲線的指標
    fpr, tpr, thresholds = roc_curve(y_test, y_proba)

    # 計算 AUC 分數
    roc_auc = roc_auc_score(y_test, y_proba)
    print(f"ROC AUC Score: {roc_auc}")

    # 繪製 ROC 曲線
    plt.figure(figsize=(8, 8))
    plt.plot(fpr, tpr, label=f"AUC = {roc_auc:.2f}")
    plt.plot([0, 1], [0, 1], "k--")  # 對角線
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel("False Positive Rate (FPR)")
    plt.ylabel("True Positive Rate (TPR)")
    plt.title("Receiver Operating Characteristic (ROC) Curve")
    plt.legend(loc="lower right")
    # plt.show()


def concatenate_csv_files(file_paths):
    """
    合併多個 CSV 檔案成一個 DataFrame。

    Parameters:
        file_paths (list): 包含多個 CSV 檔案路徑的列表。

    Returns:
        pd.DataFrame: 合併後的 DataFrame。
    """
    # 初始化一個空的 DataFrame
    concatenated_df = pd.DataFrame()

    # 迴圈讀取每個檔案並使用 pd.concat 整合到 DataFrame
    for file_path in file_paths:
        # 使用 pd.read_csv 讀取每個檔案
        current_df = pd.read_csv(file_path)

        # 使用 pd.concat 將目前的 DataFrame 與之前的合併
        concatenated_df = pd.concat([concatenated_df, current_df], ignore_index=True)

    return concatenated_df


def concatenate_csv_files_dask(file_paths):
    """
    合併多個 CSV 檔案成一個 Dask DataFrame。

    Parameters:
        file_paths (list): 包含多個 CSV 檔案路徑的列表。

    Returns:
        dask.dataframe.DataFrame: 合併後的 Dask DataFrame。
    """
    # 初始化一個空的 Dask DataFrame
    dask_dfs = []

    # 迴圈讀取每個檔案並將數據添加到 Dask DataFrame 列表
    for file_path in file_paths:
        # 使用 dask.dataframe.read_csv 讀取每個檔案，返回一個 Dask DataFrame
        current_dask_df = dd.read_csv(file_path)
        dask_dfs.append(current_dask_df)

    # 使用 dask.dataframe.concat 將 Dask DataFrame 列表合併成一個大 Dask DataFrame
    concatenated_dask_df = dd.concat(dask_dfs, axis=0)

    return concatenated_dask_df


def calculate_gdf_area(gdf_1, gdf_2, how="intersection", col="sample"):
    # Overlay the two GeoDataFrames
    temp = gpd.overlay(gdf_1, gdf_2, how=how)

    # Calculate road area grouped by village
    temp_area = temp.groupby("village")["geometry"].apply(lambda x: x.unary_union.area)
    gdf_1[col] = gdf_1["village"].map(temp_area)
    gdf_1[col] = gdf_1[col].fillna(0)
    return gdf_1


def dissolve_gdf(gdf, keep_col=["id", "geometry"]):
    """
    Dissolve a GeoDataFrame based on geometry.

    Parameters:
    - gdf (GeoDataFrame): The input GeoDataFrame to be dissolved.
    - keep_col (list): A list of column names to keep in the dissolved GeoDataFrame.

    Returns:
    - GeoDataFrame: The dissolved GeoDataFrame.
    """
    # Check missing values
    gdf = gdf[keep_col].copy()
    print(gdf.isnull().sum())

    # Check polygon validity and fix
    gdf["geometry"] = gdf["geometry"].buffer(0)

    # Dissolve
    dissolved_gdf = gdf.dissolve()
    print(dissolved_gdf.head())
    return dissolved_gdf


def calculate_pavement_area(gdf_1, gdf_2, column="area"):
    """
    Calculate the area of intersections between polygons in one GeoDataFrame and another.

    Parameters:
    gdf_1 : GeoDataFrame
        The GeoDataFrame containing the polygons.
    gdf_2 : GeoDataFrame
        The GeoDataFrame containing other geometries to calculate the intersection area with the polygons.
    column : str, default 'area'
        The name of the new column to store the intersection area. Default is 'area'.
    Returns:
    GeoDataFrame
        A GeoDataFrame containing the original polygon data along with a new column for the intersection area.
    """
    # Check if 'id' column exists before dropping it
    if "id" in gdf_2.columns:
        gdf_2 = gdf_2.rename(columns={"id": "id_right"})

    # Set crs
    gdf_1 = gdf_1.to_crs(epsg=3826)
    gdf_2 = gdf_2.to_crs(epsg=3826)

    # Calculate gdf_2 area
    gdf_2_area = gdf_2.copy()
    gdf_2_area[column] = gdf_2_area["geometry"].area

    # Spatial join
    keep_col = [column, "geometry"]
    drop_col = [col for col in gdf_2_area.columns.tolist() if col not in keep_col]
    temp = gpd.sjoin(gdf_1, gdf_2_area, how="left", predicate="intersects")
    temp = temp.drop(columns=drop_col)
    temp = temp.drop(columns=["index_right"])
    temp[column] = temp[column].fillna(0)
    temp = temp.drop_duplicates(subset=["id"])
    temp = temp.reset_index(drop=True)

    return temp


def under_sampling_negative_by_ratio(X_train, y_train, ratio=200000, random_state=42):
    """
    透過比例下採樣負樣本。

    Parameters:
    - ratio: 負樣本數量與正樣本數量的比例。
    - X_train: 訓練集的特徵。
    - y_train: 訓練集的目標變數。

    Returns:
    - X_train_n: 下採樣後的訓練集特徵。
    - y_train_n: 下採樣後的訓練集目標變數。
    """
    # Set negative over positive ratio
    n_p_ratio = ratio

    # Reset X index (only X have duplicated index)
    X_train = X_train.reset_index(drop=True)
    y_train = y_train.reset_index(drop=True)

    # Find positive and negative counts
    p_train = y_train.value_counts()[1]
    n_train = y_train.value_counts()[0]

    # List negative index
    n_train_id = y_train[y_train == 0].index
    n_train_id_df = n_train_id.to_frame()

    # Random select negative index
    num = n_train - (n_p_ratio * p_train)
    n_train_id_sampled = n_train_id_df.sample(n=num, random_state=random_state)

    # transform dataframe back to index
    n_train_id_sampled_index = n_train_id_sampled.index

    # Drop selected index on X_train and y_train
    X_train_n = X_train.drop(n_train_id_sampled_index)
    y_train_n = y_train.drop(n_train_id_sampled_index)

    # Check negative over positive ratio
    p_train_new = y_train_n.value_counts()[1]
    n_train_new = y_train_n.value_counts()[0]
    n_p_ratio_new = n_train_new / p_train_new
    print(
        f"Positive: {p_train_new}, Negative: {n_train_new}",
        "\nNegative over Positive Ratio:",
        n_p_ratio_new,
    )

    return X_train_n, y_train_n
