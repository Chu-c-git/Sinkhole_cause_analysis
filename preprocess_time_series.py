"""
This script preprocesses various environmental data for the analysis of sinkhole causes. The data 
includes tide levels, road cases, earthquake occurrences, precipitation, river levels, and 
groundwater levels. The script sets up paths and configurations, extracts and preprocesses the data,
creates time series tables, concatenates them into a single DataFrame, and exports the final 
merged DataFrame.

# Configuration
START_DATE : 
The start date for the data preprocessing period. Data before this date will be excluded.
END_DATE : 
The end date for the data preprocessing period. Data after this date will be excluded.
TIME_PERIOD : 
The time period for aggregating the data. 'M' stands for monthly aggregation.
YEARS :  
The years for which the data will be extracted and processed.
"""

import os
import xml.etree.ElementTree as ET
from datetime import datetime

import geopandas as gpd
from function import *

# Set path
DOWNLOAD_FOLDER = "data/raw_data"
EXPORT_FOLDER = "data/preprocessed_data"

# Config
START_DATE = "2019-01-01"
END_DATE = "2022-12-31"
TIME_PERIOD = "M"
YEARS = [2018, 2019, 2020, 2021, 2022, 2023]

# Function
## Data Extraction
### Tide
def extract_tide_data(download_folder):
    """
    Extract tide data from a JSON file.

    Parameters:
    download_folder (str): Path to the folder containing the JSON file.

    Returns:
    pd.DataFrame: DataFrame containing the extracted tide data.
    """
    tide_path = os.path.join(download_folder, "C-B0052-001.json")
    tide_raw = pd.read_json(tide_path, encoding="utf-8")

    # Data Overview
    station_data = tide_raw["cwaopendata"]["Resources"]["Resource"]["Data"][
        "SeaSurfaceObs"
    ]["Location"]

    # Data Extraction
    station_num = len(station_data)
    tide_record = []

    for station in range(1, station_num):
        data_year_list = station_data[station]["StationObsStatistics"]["DataYear"]
        month_num = len(station_data[station]["StationObsStatistics"]["Monthly"])

        for month in range(12, month_num):
            year_index = (month // 12) - 1
            temp = {
                "StationID": station_data[station]["Station"]["StationID"],
                "StationName": station_data[station]["Station"]["StationName"],
                "StationNameEN": station_data[station]["Station"]["StationNameEN"],
                "StationLatitude": station_data[station]["Station"]["StationLatitude"],
                "StationLongitude": station_data[station]["Station"][
                    "StationLongitude"
                ],
                "StationAttribute": station_data[station]["Station"][
                    "StationAttribute"
                ],
                "Description": station_data[station]["Station"]["Description"],
                "County": station_data[station]["Station"]["County"]["CountyName"],
                "Town": station_data[station]["Station"]["Town"]["TownName"],
                "DataYear": data_year_list[year_index],
                "DataMonth": station_data[station]["StationObsStatistics"]["Monthly"][
                    month
                ]["DataMonth"],
                "HighestHighWaterLevel": station_data[station]["StationObsStatistics"][
                    "Monthly"
                ][month]["HighestHighWaterLevel"],
                "HighestAstronomicalTide": station_data[station][
                    "StationObsStatistics"
                ]["Monthly"][month]["HighestAstronomicalTide"],
                "MeanHighWaterLevel": station_data[station]["StationObsStatistics"][
                    "Monthly"
                ][month]["MeanHighWaterLevel"],
                "MeanTideLevel": station_data[station]["StationObsStatistics"][
                    "Monthly"
                ][month]["MeanTideLevel"],
                "MeanLowWaterLevel": station_data[station]["StationObsStatistics"][
                    "Monthly"
                ][month]["MeanLowWaterLevel"],
                "LowestAstronomicalTide": station_data[station]["StationObsStatistics"][
                    "Monthly"
                ][month]["LowestAstronomicalTide"],
                "LowestLowWaterLevel": station_data[station]["StationObsStatistics"][
                    "Monthly"
                ][month]["LowestLowWaterLevel"],
                "MeanTidalRange": station_data[station]["StationObsStatistics"][
                    "Monthly"
                ][month]["MeanTidalRange"],
                "MaxAstronomicalTidalRange": station_data[station][
                    "StationObsStatistics"
                ]["Monthly"][month]["MaxAstronomicalTidalRange"],
                "MeanHighWaterOfSpringTide": station_data[station][
                    "StationObsStatistics"
                ]["Monthly"][month]["MeanHighWaterOfSpringTide"],
                "MeanLowWaterOfSpringTide": station_data[station][
                    "StationObsStatistics"
                ]["Monthly"][month]["MeanLowWaterOfSpringTide"],
            }
            tide_record.append(temp)

    tide_df = pd.DataFrame(tide_record)
    return tide_df


### Road case
def extract_road_case_data(download_folder):
    """
    Extract road case data from a CSV file.

    Parameters:
    download_folder (str): Path to the folder containing the CSV file.

    Returns:
    pd.DataFrame: DataFrame containing the extracted road case data.
    """
    road_case_path = os.path.join(
        download_folder, "道管系統坑洞案件_108-112_Chu加案件標記_20240201.csv"
    )

    # Check if the file exists
    if os.path.exists(road_case_path):

        print("Road case data found, proceed with further data preprocessing.")
        road_case_raw = pd.read_csv(road_case_path, encoding="utf-8")

        # Data Preprocessing
        road_case_raw["查報日期"] = pd.to_datetime(road_case_raw["查報日期"])

        # Filter columns
        keep_col = ["案件編號", "查報日期", "lng_twd97", "lat_twd97", "TUIC天坑判斷"]
        road_case_df = road_case_raw[keep_col]

        # Mask
        mask_sinkhole = road_case_df["TUIC天坑判斷"] == 2
        road_case_df = road_case_df[mask_sinkhole]

        return road_case_df
    else:
        print(
            "No road case data found, skip and proceed with other time series data processing."
        )
        return None


### Taipei border
def extract_tp_border_data(download_folder):
    """
    Extract Taipei border data from a shapefile.

    Parameters:
    download_folder (str): Path to the folder containing the shapefile.

    Returns:
    gpd.GeoDataFrame: GeoDataFrame containing the extracted Taipei border data.
    """
    tp_border_path = os.path.join(
        download_folder, "臺北市區界圖_20220915", "G97_A_CADIST_P.shp"
    )
    tp_border_raw = gpd.read_file(tp_border_path, encoding="utf-8")

    # Set tp crs
    tp_border_raw.set_crs(epsg=3826, inplace=True)

    # Dissolve
    tp_border_gdf = tp_border_raw.dissolve()
    tp_border_gdf["geometry"] = tp_border_gdf["geometry"].buffer(distance=100 * 1000)

    return tp_border_gdf


### Earthquake data
def extract_earthquake_data(download_folder, tp_border_gdf, years):
    """
    Extract earthquake data from a XML file.

    Parameters:
    download_folder (str): Path to the folder containing the XML file.

    Returns:
    pd.DataFrame: DataFrame containing the extracted earthquake data.
    """
    earthquake_df = pd.DataFrame()

    for year in years:
        xml_path = os.path.join(
            download_folder, "E-A0073-002", f"CWA-EQ-Catalog-{year}.xml"
        )
        # Read XML
        tree = ET.parse(xml_path)
        root = tree.getroot()
        df_list = []
        for i in range(2, len(root[8][1])):
            data = root[8][1][i]
            if len(data) == 14:
                temp = {
                    "發震時間(OriginTime)": data[0].text,
                    "震央經度(EpicenterLng)": float(data[1].text),
                    "震央緯度(EpicenterLat)": float(data[2].text),
                    "震源深度(FocalDepth)": data[3].text,
                    "芮氏規模(LocalMagnitude)": data[4].text,
                    "定位測站個數(StationNumber)": data[5].text,
                    "定位相位個數(PhaseNumber)": data[6].text,
                    "震央距(MinimumEpicenterDistance)": data[7].text,
                    "最大空餘角度(gap)": data[8].text,
                    "震波走時殘差(rms)": data[9].text,
                    "水平標準偏差(erh)": data[10].text,
                    "垂直標準偏差(erz)": data[11].text if len(data) > 13 else None,
                    "定位品質(Quality)": (
                        data[12].text if len(data) > 13 else data[11].text
                    ),
                    "定位回顧狀態(ReviewStatus)": (
                        data[13].text if len(data) > 13 else data[12].text
                    ),
                }
                df_list.append(temp)
            elif len(data) == 13:
                temp = {
                    "發震時間(OriginTime)": data[0].text,
                    "震央經度(EpicenterLng)": float(data[1].text),
                    "震央緯度(EpicenterLat)": float(data[2].text),
                    "震源深度(FocalDepth)": data[3].text,
                    "芮氏規模(LocalMagnitude)": data[4].text,
                    "定位測站個數(StationNumber)": data[5].text,
                    "定位相位個數(PhaseNumber)": data[6].text,
                    "震央距(MinimumEpicenterDistance)": data[7].text,
                    "最大空餘角度(gap)": data[8].text,
                    "震波走時殘差(rms)": data[9].text,
                    "水平標準偏差(erh)": data[10].text,
                    "垂直標準偏差(erz)": None,
                    "定位品質(Quality)": data[11].text,
                    "定位回顧狀態(ReviewStatus)": data[12].text,
                }
                df_list.append(temp)
            elif len(data) == 12:
                temp = {
                    "發震時間(OriginTime)": data[0].text,
                    "震央經度(EpicenterLng)": float(data[1].text),
                    "震央緯度(EpicenterLat)": float(data[2].text),
                    "震源深度(FocalDepth)": data[3].text,
                    "芮氏規模(LocalMagnitude)": data[4].text,
                    "定位測站個數(StationNumber)": data[5].text,
                    "定位相位個數(PhaseNumber)": data[6].text,
                    "震央距(MinimumEpicenterDistance)": data[7].text,
                    "最大空餘角度(gap)": data[8].text,
                    "震波走時殘差(rms)": data[9].text,
                    "水平標準偏差(erh)": None,
                    "垂直標準偏差(erz)": None,
                    "定位品質(Quality)": data[10].text,
                    "定位回顧狀態(ReviewStatus)": data[11].text,
                }
                df_list.append(temp)
            elif len(data) < 12:
                temp = {
                    "發震時間(OriginTime)": data[0].text,
                    "震央經度(EpicenterLng)": float(data[1].text),
                    "震央緯度(EpicenterLat)": float(data[2].text),
                    "震源深度(FocalDepth)": data[3].text,
                    "芮氏規模(LocalMagnitude)": data[4].text,
                    "定位測站個數(StationNumber)": data[5].text,
                    "定位相位個數(PhaseNumber)": data[6].text,
                    "震央距(MinimumEpicenterDistance)": data[7].text,
                    "最大空餘角度(gap)": data[8].text,
                    "震波走時殘差(rms)": None,
                    "水平標準偏差(erh)": None,
                    "垂直標準偏差(erz)": None,
                    "定位品質(Quality)": data[9].text,
                    "定位回顧狀態(ReviewStatus)": data[10].text,
                }
                df_list.append(temp)
            df_year = pd.DataFrame(df_list)
        earthquake_df = pd.concat([earthquake_df, df_year], axis=0)
        # Spatial join
        gdf_earthquake = gpd.GeoDataFrame(
            earthquake_df,
            geometry=gpd.points_from_xy(
                earthquake_df["震央經度(EpicenterLng)"],
                earthquake_df["震央緯度(EpicenterLat)"],
            ),
        )
        gdf_earthquake = gdf_earthquake.set_crs(epsg=4326)
        gdf_earthquake = gdf_earthquake.to_crs(epsg=3826)
        earthquake_gdf = gpd.sjoin(
            gdf_earthquake, tp_border_gdf, how="inner", predicate="intersects"
        )
        # Drop columns
        drop_col = tp_border_gdf.columns.to_list()
        earthquake_gdf = earthquake_gdf.copy()
        earthquake_gdf.drop(columns=drop_col, inplace=True)

    return earthquake_gdf



# Extract all rainfall data
def extract_all_rainfall_data(download_folder, years):
    """
    Extract monthly rainfall data for each station from multiple XML files and return a DataFrame.

    Parameters:
    years (list): A list of years for which the data should be extracted.
    download_folder (str): Path to the folder containing the XML files.

    Returns:
    pd.DataFrame: DataFrame containing the monthly rainfall data for each station for all years.
    """

    ### Monthly precipitation of a year
    def extract_yearly_rainfall_data(root, num_location_nodes):
        """
        Extract monthly rainfall data for each station from an XML tree and return a DataFrame.

        Parameters:
        root (Element): The root node of the XML tree.
        num_location_nodes (int): The number of location nodes.

        Returns:
        pd.DataFrame: DataFrame containing the monthly rainfall data for each station.
        """
        # Initialize an empty list to store all the data
        data = []

        # Iterate over each station and collect the data
        for sta in range(num_location_nodes):
            ch_name = root[10][0][1][0][sta][0][1].text
            station_id = root[10][0][1][0][sta][0][0].text
            en_name = root[10][0][1][0][sta][0][2].text
            attr = root[10][0][1][0][sta][0][3].text

            # Put the station data into the list, including the data for each month
            data.extend(
                [
                    {
                        "ch_name": ch_name,
                        "station_id": station_id,
                        "en_name": en_name,
                        "attr": attr,
                        "year_month": root[10][0][1][0][sta][2][0][b][0].text,
                        "total_rain": root[10][0][1][0][sta][2][0][b][1].text,
                    }
                    for b in range(12)
                ]
            )

        # Generate a DataFrame using the collected data
        cwb_data = pd.DataFrame(data)
        return cwb_data
    
    # Initialize an empty list to store all the data
    all_data = []

    # Iterate over each year and extract the data
    for year in years:
        xml_path = os.path.join(download_folder, "C-B0025-002", f"dy_Report_{year}.xml")

        # Read XML
        tree = ET.parse(xml_path)
        root = tree.getroot()

        # Count the number of location nodes
        namespace = {"ns": "urn:cwa:gov:tw:cwacommon:0.1"}
        location_nodes = root.findall(".//ns:location", namespace)
        num_location_nodes = len(location_nodes)

        # Extract the data for the current year
        year_data = extract_yearly_rainfall_data(root, num_location_nodes)

        # Append the data to the list
        all_data.append(year_data)

    # Concatenate all the data into a single DataFrame
    precipitation_df = pd.concat(all_data, ignore_index=True)
    return precipitation_df


### Extract river level data
def extract_river_level_data(download_folder):
    """
    Extract river level data from a CSV file.

    Parameters:
    download_folder (str): Path to the folder containing the CSV file.

    Returns:
    pd.DataFrame: DataFrame containing the extracted river level data.
    """
    river_level_path = os.path.join(download_folder, "202400337_日水位.csv")
    river_level_raw = pd.read_csv(river_level_path, encoding="big5")

    # Transform wide sheet to long sheet
    river_data = pd.melt(
        river_level_raw,
        id_vars=["管理單位", "站名", "站號", "年份", "月份"],
        var_name="日期",
        value_name="水位(m)",
    )

    # Clean '日期' column
    river_data["日期"] = river_data["日期"].str.replace("日(m)", "")

    # Add Date column
    str_col = ["年份", "月份", "日期"]
    river_data[str_col] = river_data[str_col].astype(str)
    river_data["date"] = (
        river_data["年份"] + "-" + river_data["月份"] + "-" + river_data["日期"]
    )

    # Transform date and sortby date
    river_data["date"] = pd.to_datetime(
        river_data["date"], format="%Y-%m-%d", errors="coerce"
    )
    river_data["date"] = river_data["date"].dt.date
    river_data = river_data.sort_values(by="date")

    river_data.dropna(subset=["date"], inplace=True)
    river_data = river_data.reset_index(drop=True)
    
    # river level station
    river_station_path = os.path.join(
        DOWNLOAD_FOLDER, "RIVWLSTA_e", "RIVWLSTA_e", "RIVWLSTA_e.shp"
    )
    river_station_raw = gpd.read_file(river_station_path, encoding="utf-8")

    # Set crs
    river_station_raw = river_station_raw.set_crs(epsg=3826, inplace=True)
    river_station_raw = river_station_raw.to_crs(epsg=4326)

    # Filter columns
    keep_col = ["ST_NO", "geometry"]
    river_station = river_station_raw[keep_col]
    river_station.head()
    # Strip columns
    river_data["站號"] = river_data["站號"].str.strip()

    # Merge river data with river station
    river_data = river_data.merge(
        river_station, left_on="站號", right_on="ST_NO", how="left"
    )

    # Drop columns
    river_data = river_data.drop(columns=["ST_NO"])

    # Fill missing value
    river_data["水位(m)"] = river_data["水位(m)"].replace("缺測", None)

    # Transform geometry
    river_level_gdf = gpd.GeoDataFrame(river_data, geometry="geometry", crs="EPSG:4326")

    # Sort value by station and date
    river_level_gdf = river_level_gdf.sort_values(by=["站名", "date"], ascending=True)
    river_level_gdf.reset_index(inplace=True)
    river_level_gdf = river_level_gdf.drop(columns="index")

    return river_level_gdf


### Extract Ground water level data
def extract_groundwater_level_data(download_folder):
    """
    Extract ground water level data from a CSV file.

    Parameters:
    download_folder (str): Path to the folder containing the CSV file.

    Returns:
    pd.DataFrame: DataFrame containing the extracted ground water level data.
    """
    groundwater_level_path = os.path.join(download_folder, "202400337_自記站日水位.csv")
    groundwater_level_raw = pd.read_csv(groundwater_level_path, encoding="big5")

    groundwater_level_data = pd.melt(
        groundwater_level_raw,
        id_vars=["管理單位", "井名", "井號", "年份", "月份"],
        var_name="日期",
        value_name="水位(m)",
    )

    # Clean '日期' column
    groundwater_level_data["日期"] = groundwater_level_data["日期"].str.replace("日(m)", "")

    # Add Date column
    str_col = ["年份", "月份", "日期"]
    groundwater_level_data[str_col] = groundwater_level_data[str_col].astype(str)
    groundwater_level_data["date"] = (
        groundwater_level_data["年份"]
        + "-"
        + groundwater_level_data["月份"]
        + "-"
        + groundwater_level_data["日期"]
    )

    # Transform date and sortby date
    groundwater_level_data["date"] = pd.to_datetime(
        groundwater_level_data["date"], format="%Y-%m-%d", errors="coerce"
    )
    groundwater_level_data["date"] = groundwater_level_data["date"].dt.date
    groundwater_level_data = groundwater_level_data.sort_values(by="date")
    groundwater_level_data.dropna(subset=["date"], inplace=True)
    groundwater_level_data = groundwater_level_data.reset_index(drop=True)

    # Fill missing value
    groundwater_level_data["水位(m)"] = groundwater_level_data["水位(m)"].replace("缺測", None)
    groundwater_level_data["水位(m)"] = (
        groundwater_level_data["水位(m)"].astype(float).round(2)
    )

    # Ground water station
    groundwater_station_path = os.path.join(
        DOWNLOAD_FOLDER, "gwobwell_e", "gwobwell_e", "gwobwell_e.shp"
    )
    groundwater_station_raw = gpd.read_file(groundwater_station_path, encoding="utf-8")

    # Filter columns
    keep_col = ["井名", "geometry"]
    groundwater_station = groundwater_station_raw[keep_col]

    # Transform geometry
    groundwater_station.loc[:, "geometry"] = groundwater_station.loc[:, "geometry"].set_crs(
        epsg=3826
    )
    groundwater_station.loc[:, "geometry"] = groundwater_station.loc[:, "geometry"].to_crs(
        epsg=4326
    )

    # Find geometry from station name
    groundwater_level_data["井名"] = groundwater_level_data["井名"].str.strip()
    groundwater_level_data = groundwater_level_data.merge(
        groundwater_station, how="left", left_on="井名", right_on="井名"
    )

    # Transform df to gdf
    groundwater_level_gdf = gpd.GeoDataFrame(groundwater_level_data, geometry="geometry")

    # Sort value by well and date
    groundwater_level_gdf = groundwater_level_gdf.sort_values(
        by=["井名", "date"], ascending=True
    )
    groundwater_level_gdf.reset_index(inplace=True)
    groundwater_level_gdf = groundwater_level_gdf.drop(columns="index")

    return groundwater_level_gdf


# Build Time Series Table
## Tide level
def create_tide_table(tide_df, start_date, end_date, time_period):
    """
    Create a time series table for tide level data.

    Parameters:
    tide_df (pd.DataFrame): DataFrame containing the tide level data.
    start_date (str): Start date for the time series table.
    end_date (str): End date for the time series table.
    time_period (str): Time period for the time series table.

    Returns:
    pd.DataFrame: Time series table for tide level data.
    """
    # Create date column
    tide_df["DataYear"] = tide_df["DataYear"].astype(int)
    tide_df["DataMonth"] = tide_df["DataMonth"].astype(int)
    tide_df["Date"] = tide_df.apply(
        lambda row: datetime(row["DataYear"], row["DataMonth"], 1), axis=1
    )

    # Filter columns
    keep_col_t = [
        "Date",
        "StationName",
        "MeanTideLevel",
        "MeanHighWaterLevel",
        "MeanLowWaterLevel",
    ]
    tide_cl = tide_df[keep_col_t]

    # mask
    mask_date = (tide_cl["Date"] >= start_date) & (tide_cl["Date"] <= end_date)
    mask_station = tide_cl["StationName"] == "淡水潮位站"
    mask_tdie = mask_date & mask_station
    tide_cl = tide_cl[mask_tdie]

    # Filter columns
    period_count = tide_cl.copy()
    period_count["Date"] = period_count["Date"].dt.to_period(time_period)

    # period_count_tide = period_count[['Date', 'MeanTideLevel']].reset_index(drop=True)
    period_count_tide = period_count.drop(columns=["StationName"])

    period_count_tide = period_count_tide.set_index("Date")
    return period_count_tide

## Road case
def create_road_case_table(road_case_df, start_date, end_date, time_period):
    """
    Create a time series table for road case data.

    Parameters:
    road_case_df (pd.DataFrame): DataFrame containing the road case data.
    start_date (str): Start date for the time series table.
    end_date (str): End date for the time series table.
    time_period (str): Time period for the time series table.

    Returns:
    pd.DataFrame: Time series table for road case data.
    """
    if road_case_df is not None:
        # Filter columns
        keep_col_r = ["查報日期", "lng_twd97", "lat_twd97"]
        road_case_cl = road_case_df[keep_col_r]

        # mask
        mask_date = (road_case_cl["查報日期"] >= start_date) & (
            road_case_cl["查報日期"] <= end_date
        )
        mask_sinkhole = road_case_df["TUIC天坑判斷"] == 2
        mask_road_case = mask_date & mask_sinkhole
        road_case_cl = road_case_cl[mask_road_case]

        # Roadcase groupby
        period_count = road_case_cl.copy()
        period_count["查報日期"] = period_count["查報日期"].dt.to_period(time_period)
        period_count_sinkhole = period_count.groupby("查報日期").size().reset_index()
        period_count_sinkhole = period_count_sinkhole.rename(columns={0: "sinkhole_count"})
        period_count_sinkhole = period_count_sinkhole.set_index("查報日期")
        return period_count_sinkhole

    else:
        period_count_sinkhole = None
        print("Skipping road case data processing.")


## Earthquake
def create_earthquake_table(earthquake_gdf, start_date, end_date, time_period):
    """
    Create a time series table for earthquake data.

    Parameters:
    earthquake_gdf (gpd.GeoDataFrame): GeoDataFrame containing the earthquake data.
    start_date (str): Start date for the time series table.
    end_date (str): End date for the time series table.
    time_period (str): Time period for the time series table.

    Returns:
    pd.DataFrame: Time series table for earthquake data.
    """
    # Date
    earthquake_gdf["發震時間(OriginTime)"] = pd.to_datetime(
        earthquake_gdf["發震時間(OriginTime)"]
    )
    earthquake_gdf["芮氏規模(LocalMagnitude)"] = earthquake_gdf[
        "芮氏規模(LocalMagnitude)"
    ].astype(float)
    earthquake_gdf["震源深度(FocalDepth)"] = earthquake_gdf["震源深度(FocalDepth)"].astype(
        float
    )

    # Filter columns
    keep_col_e = [
        "發震時間(OriginTime)",
        "芮氏規模(LocalMagnitude)",
        "震源深度(FocalDepth)",
    ]
    earthquake_cl = earthquake_gdf[keep_col_e]

    # Mask
    mask_date = (earthquake_cl["發震時間(OriginTime)"] >= start_date) & (
        earthquake_cl["發震時間(OriginTime)"] <= end_date
    )
    mask_mag = earthquake_cl["芮氏規模(LocalMagnitude)"] >= 5
    mask_depth = earthquake_cl["震源深度(FocalDepth)"] <= 70
    mask_earthquake = mask_date & mask_mag & mask_depth
    earthquake_cl = earthquake_cl[mask_earthquake]

    # Resample & count
    period_count = earthquake_cl.copy()
    period_count["發震時間(OriginTime)"] = period_count[
        "發震時間(OriginTime)"
    ].dt.tz_localize(None).dt.to_period(time_period)
    period_count_earthquake = (
        period_count.groupby("發震時間(OriginTime)").size().reset_index()
    )
    period_count_earthquake = period_count_earthquake.rename(
        columns={0: "earthquake_count"}
    )
    period_count_earthquake = period_count_earthquake.set_index("發震時間(OriginTime)")
    return period_count_earthquake

## period precipitation
def create_precipitation_table(precipitation_df, start_date, end_date, time_period):
    """
    Create a time series table for precipitation data.

    Parameters:
    precipitation_df (pd.DataFrame): DataFrame containing the precipitation data.
    start_date (str): Start date for the time series table.
    end_date (str): End date for the time series table.
    time_period (str): Time period for the time series table.

    Returns:
    pd.DataFrame: Time series table for precipitation data.
    """
    # Date
    precipitation_df["year_month"] = pd.to_datetime(precipitation_df["year_month"])

    # Keep columns
    keep_col_pr = ["ch_name", "year_month", "total_rain"]
    precipitation_cl = precipitation_df[keep_col_pr]

    # mask
    mask_date = (precipitation_cl["year_month"] >= start_date) & (
        precipitation_cl["year_month"] <= end_date
    )
    mask_station = precipitation_cl["ch_name"] == "臺北"
    mask_precipitation = mask_date & mask_station
    precipitation_cl = precipitation_cl[mask_precipitation]
    precipitation_cl = precipitation_cl.drop(columns=["ch_name"])

    # Filter columns
    period_count = precipitation_cl.copy()
    period_count["year_month"] = period_count["year_month"].dt.to_period(time_period)
    period_count_precipitation = period_count.set_index("year_month")
    return period_count_precipitation

## River level
def create_river_level_table(river_level_gdf, start_date, end_date, time_period):
    """
    Create a time series table for river level data.

    Parameters:
    river_level_gdf (gpd.GeoDataFrame): GeoDataFrame containing the river level data.
    start_date (str): Start date for the time series table.
    end_date (str): End date for the time series table.
    time_period (str): Time period for the time series table.

    Returns:
    pd.DataFrame: Time series table for river level data.
    """
    # Fill missing value
    river_level_gdf["水位(m)"] = river_level_gdf["水位(m)"].astype(float)
    river_level_gdf = fill_missing_level(river_level_gdf)
    river_level_gdf["水位(m)"] = river_level_gdf["水位(m)"].round(2)

    # convert to datetime
    river_level_gdf["date"] = pd.to_datetime(river_level_gdf["date"])

    # Strip
    river_level_gdf["站名"] = river_level_gdf["站名"].str.strip()

    # Keep columns
    keep_col_river = ["date", "站名", "水位(m)"]
    river_level_cl = river_level_gdf[keep_col_river]

    # mask
    mask_date = (river_level_gdf["date"] >= start_date) & (
        river_level_gdf["date"] <= end_date
    )
    river_level_cl = river_level_cl[mask_date]

    river_level_cl.head()
    # Change to pivot table
    df = river_level_cl.copy()
    temp_river_level = df.pivot_table(index=["date"], columns="站名", values="站名")
    temp_river_level.reset_index(inplace=True)

    # Rename columns
    col_dict = {
        "寶橋": "river_level_baoqiao",
        "萬福橋": "river_level_wanfu",
    }
    temp_river_level.rename(columns=col_dict, inplace=True)
    temp_river_level.index.name = None
    temp_river_level.head()
    # Convert to datetime
    temp_river_level["date"] = pd.to_datetime(temp_river_level["date"])

    # Set index
    temp_river_level.set_index("date", inplace=True)
    period_river_level = temp_river_level.resample(time_period).median()
    period_river_level.index = period_river_level.index.to_period(time_period)
    period_river_level = period_river_level.round(2)

    # Calculate mean
    period_river_level["river_level_mean"] = period_river_level.mean(axis=1)
    period_river_level = period_river_level.round(2)

    return period_river_level


## Ground water level
def create_groundwater_level_table(
    groundwater_level_gdf, start_date, end_date
    ):
    """
    Create a time series table for ground water level data.

    Parameters:
    groundwater_level_gdf (gpd.GeoDataFrame): GeoDataFrame containing the ground water level data.
    start_date (str): Start date for the time series table.
    end_date (str): End date for the time series table.
    time_period (str): Time period for the time series table.

    Returns:
    pd.DataFrame: Time series table for ground water level data.
    """
    # Fill missing level
    groundwater_level_gdf["水位(m)"] = groundwater_level_gdf["水位(m)"].astype(float)
    groundwater_level_gdf = fill_missing_level(groundwater_level_gdf)
    groundwater_level_gdf["水位(m)"] = groundwater_level_gdf["水位(m)"].round(2)

    # Convert to datetime
    groundwater_level_gdf["date"] = pd.to_datetime(groundwater_level_gdf["date"])

    # Strip
    groundwater_level_gdf["井名"] = groundwater_level_gdf["井名"].str.strip()

    # Keep columns
    keep_col_ugwater = ["date", "井名", "水位(m)"]
    ugwater_level_cl = groundwater_level_gdf[keep_col_ugwater]

    # mask
    mask_date = (groundwater_level_gdf["date"] >= start_date) & (
        groundwater_level_gdf["date"] <= end_date
    )
    ugwater_level_cl = ugwater_level_cl[mask_date]
    ugwater_level_cl.head()
    # Change to pivot table
    df = ugwater_level_cl.copy()
    temp_ugwater_level = df.pivot_table(index=["date"], columns="井名", values="井名")
    temp_ugwater_level.reset_index(inplace=True)

    # Rename columns
    col_dict = {
        "北投(1)": "ugwater_level_bt1",
        "北投(2)": "ugwater_level_bt2",
        "台大(1)": "ugwater_level_ntu1",
        "台大(2)": "ugwater_level_ntu2",
        "國父紀念館": "ugwater_level_sunn",
        "大橋國小": "ugwater_level_dq",
        "敦化南路": "ugwater_level_dh",
        "新生公園": "ugwater_level_xs",
        "青年公園": "ugwater_level_qn",
    }
    temp_ugwater_level.rename(columns=col_dict, inplace=True)
    temp_ugwater_level.head()
    # Convert to datetime
    temp_ugwater_level["date"] = pd.to_datetime(temp_ugwater_level["date"])

    # Set index
    temp_ugwater_level.set_index("date", inplace=True)
    period_ugwater_level = temp_ugwater_level.resample(TIME_PERIOD).median()
    period_ugwater_level.index = period_ugwater_level.index.to_period(TIME_PERIOD)
    period_ugwater_level = period_ugwater_level.round(2)

    # Calculate mean
    period_ugwater_level["ugwater_level_mean"] = period_ugwater_level.mean(axis=1)
    period_ugwater_level = period_ugwater_level.round(2)
    
    return period_ugwater_level

# Concate Dataframe
## concat dataframes
def concat_dataframes(
    period_count_sinkhole,
    period_count_tide,
    period_count_earthquake,
    period_count_precipitation,
    period_river_level,
    period_ugwater_level,
):
    """
    Concatenate multiple DataFrames into a single DataFrame.

    Parameters:
    period_count_sinkhole (pd.DataFrame): DataFrame containing the sinkhole count data.
    period_count_tide (pd.DataFrame): DataFrame containing the tide level data.
    period_count_earthquake (pd.DataFrame): DataFrame containing the earthquake count data.
    period_count_precipitation (pd.DataFrame): DataFrame containing the precipitation data.
    period_river_level (pd.DataFrame): DataFrame containing the river level data.
    period_ugwater_level (pd.DataFrame): DataFrame containing the ground water level data.

    Returns:
    pd.DataFrame: Concatenated DataFrame containing all the input DataFrames.
    """
    dataframes = [
        period_count_sinkhole,
        period_count_tide,
        period_count_earthquake,
        period_count_precipitation,
        period_river_level,
        period_ugwater_level,
    ]

    # Filter out None values
    not_none_dataframes = [df for df in dataframes if df is not None]

    # Concatenate non-None DataFrames
    if not_none_dataframes:
        merged_df = pd.concat(not_none_dataframes, axis=1)
    else:
        merged_df = pd.DataFrame()  # Return an empty DataFrame if all are None

    return merged_df

## Create tide range feature
def arrange_merged_df(merged_df, time_period):
    """
    Arrange the columns of the merged DataFrame.

    Parameters:
    merged_df (pd.DataFrame): DataFrame containing the merged data.
    time_period (str): Time period for the time series table.

    Returns:
    pd.DataFrame: DataFrame with the columns arranged in the desired order.
    """
    # Sort index & fillna
    merged_df = merged_df.sort_index().fillna(0)
    merged_df = merged_df.reset_index()

    # Period to timestamp
    merged_df["index"] = merged_df["index"].dt.to_timestamp(freq=time_period)

    # Create tide range features
    merged_df["MeanHighWaterLevel"] = merged_df["MeanHighWaterLevel"].astype(float)
    merged_df["MeanLowWaterLevel"] = merged_df["MeanLowWaterLevel"].astype(float)
    merged_df["tide_range"] = (
        merged_df["MeanHighWaterLevel"] - merged_df["MeanLowWaterLevel"]
    )
    merged_df["tide_range"] = merged_df["tide_range"].fillna(0)
    return merged_df

## Export
def export_merged_df(merged_df, export_folder, time_period):
    """
    Export the merged DataFrame to a CSV file.

    Parameters:
    merged_df (pd.DataFrame): DataFrame containing the merged data.
    export_folder (str): Path to the folder where the CSV file should be saved.
    time_period (str): Time period for the time series table.
    """
    # Create export folder
    def create_export_folder():
        if not os.path.exists(EXPORT_FOLDER):
            os.makedirs(EXPORT_FOLDER)
        else:
            print("Export folder exists.")
    create_export_folder()

    # Export to CSV
    export_path = os.path.join(export_folder, f"time_series_table_{time_period}.csv")
    merged_df.to_csv(export_path, index=False, encoding="utf-8")
    print(f"Time series table saved to: {export_path}")

def main():
    # Extract data
    tide_df = extract_tide_data(DOWNLOAD_FOLDER)
    road_case_df = extract_road_case_data(DOWNLOAD_FOLDER)
    tp_border_gdf = extract_tp_border_data(DOWNLOAD_FOLDER)
    earthquake_gdf = extract_earthquake_data(DOWNLOAD_FOLDER, tp_border_gdf, YEARS)
    precipitation_df = extract_all_rainfall_data(DOWNLOAD_FOLDER, YEARS)
    river_level_gdf = extract_river_level_data(DOWNLOAD_FOLDER)
    groundwater_level_gdf = extract_groundwater_level_data(DOWNLOAD_FOLDER)

    # Create time series tables
    period_count_sinkhole = create_road_case_table(
        road_case_df, START_DATE, END_DATE, TIME_PERIOD
    )
    period_count_tide = create_tide_table(
        tide_df, START_DATE, END_DATE, TIME_PERIOD
    )
    period_count_earthquake = create_earthquake_table(
        earthquake_gdf, START_DATE, END_DATE, TIME_PERIOD
    )
    period_count_precipitation = create_precipitation_table(
        precipitation_df, START_DATE, END_DATE, TIME_PERIOD
    )
    period_river_level = create_river_level_table(
        river_level_gdf, START_DATE, END_DATE, TIME_PERIOD
    )
    period_ugwater_level = create_groundwater_level_table(
        groundwater_level_gdf, START_DATE, END_DATE, TIME_PERIOD
    )

    # Concatenate DataFrames
    merged_df = concat_dataframes(
        period_count_sinkhole,
        period_count_tide,
        period_count_earthquake,
        period_count_precipitation,
        period_river_level,
        period_ugwater_level,
    )

    # Arrange the merged DataFrame
    arranged_df = arrange_merged_df(merged_df, TIME_PERIOD)

    # Export the merged DataFrame
    export_merged_df(arranged_df, EXPORT_FOLDER, TIME_PERIOD)

if __name__ == "__main__":

    START_TIME = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f'Data Preprocessing Start at {START_TIME}!\n')

    main()

    END_TIME = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f'\nData Preprocessing Finished at {END_TIME}!')