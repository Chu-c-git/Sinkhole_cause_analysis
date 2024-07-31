"""
This script is designed to preprocess and analyze hexagon-based spatial data for Taipei City,
focusing on road intersections and geographic attributes. The primary functionalities include:
1. Checking and creating necessary directories and data files.
2. Loading and processing geospatial data, including hexagon grids and road networks.
3. Intersecting hexagon grids with road data to determine the areas covered by roads.
4. Exporting processed data for further analysis or visualization.

Configuration Parameters:
--------------------------
- download_folder (str): The directory where raw data is stored.
  Default: "data/raw_data"

- export_folder (str): The directory where preprocessed data will be saved.
  Default: "data/preprocessed_data"

- temporal_folder (str): A temporary directory for intermediate data storage.
  Default: "data/temporal_data"

- select_town (list of str): The list of towns to filter the final data.
  Default: ["大同區"]

- export_tp_all (bool): Whether to export all hexagons covering Taipei City.
  Default: False

- export_tp_road_all (bool): Whether to export all hexagons intersected with roads.
  Default: True

- export_tp_road_selected (bool): Whether to export only hexagons within the selected towns.
  Default: False (Recommended to set to True if you want to filter by town)
"""

import os
import numpy as np
import glob
import psutil
import shutil
from datetime import datetime
from shapely import wkt

import geopandas as gpd
import dask_geopandas as dgpd
from function import *

class PreprocessHex:
    def __init__(
            self, 
            download_folder="data/raw_data", 
            export_folder="data/preprocessed_data",
            temporal_folder="data/temporal_data",
            select_town=["大同區"],
            export_tp_all=False,
            export_tp_road_all=True,
            export_tp_road_selected=False   
        ):
        self.download_folder = download_folder
        self.export_folder = export_folder
        self.select_town = select_town
        self.temporal_folder = temporal_folder
        self.export_tp_all = export_tp_all
        self.export_tp_road_all = export_tp_road_all
        self.export_tp_road_selected = export_tp_road_selected

    def check_temporal_folder(self)->None:
        """
        Check if temporal folder exist, if not create one.
        """
        if not os.path.exists(self.temporal_folder):
            os.makedirs(self.temporal_folder)
        else:
            print("temporal_folder check: exist")

    def check_original_hex_file(self)->None:
        """
        Check if the original hexagon file exist
        """
        hex_5m_path = os.path.join(self.download_folder, "tp_original_hex.gpkg")
        if os.path.exists(hex_5m_path):
            print("tp_original_hex.gpkg check: exist.")
        else:
            print("tp_original_hex.gpkg not found, please make sure you have created the file by QGIS!")            

    def check_memory_and_storage(self)->bool:
        # Check if the system has at least 16GB of RAM
        total_memory_gb = psutil.virtual_memory().total / (1024 ** 3)
        if total_memory_gb < 15:
            print(f"Warning: The system memory is less than 15GB (Current memory: {total_memory_gb:.2f} GB)")
        else:
            print(f"System memory is sufficient (Current memory: {total_memory_gb:.2f} GB)")
        
        # Check if there is at least 10GB of free space on the main storage drive
        free_space_gb = shutil.disk_usage("/").free / (1024 ** 3)
        if free_space_gb < 10:
            print(f"Warning: The available storage space is less than 10GB (Current available space: {free_space_gb:.2f} GB)")
        else:
            print(f"Available storage space is sufficient (Current available space: {free_space_gb:.2f} GB)")

        # Prompt the user to continue or not
        proceed = input("Do you want to proceed with the execution? (y/n): ").strip().lower()
        if proceed == 'y':
            return True
        else:
            print("Operation cancelled.")
            return False

    def load_hex(self)->gpd.GeoDataFrame:
        """
        Load 5m hexagon data and create centroid for future use.
        Approximate time: 
        - read_file: 2min
        - find_centroid: 7min
        - compute_drop_col: 8min
        """
        hex_5m_path = os.path.join(self.download_folder, "tp_original_hex.gpkg")
        hex_5m = dgpd.read_file(hex_5m_path, chunksize=100e6)
        hex_5m["centroid"] = hex_5m["geometry"].centroid
        drop_col = ["left", "right", "top", "bottom"]
        hex_5m = hex_5m.drop(columns=drop_col)
        hex_5m = hex_5m.compute()
        return hex_5m

    def split_hex(self, hex_5m)->None:
        """
        Split hexagon data into chunks and save to temporal folder
        Approximate time: 
        - split: 34min
        """
        chunk_size = 100000
        arrange = np.arange(len(hex_5m)) // chunk_size
        count = (len(hex_5m) // chunk_size) + 1
        grouped_chunks = pd.DataFrame(hex_5m).groupby(arrange)

        for i, (_, chunk) in enumerate(grouped_chunks):
            chunk.to_csv(f'{self.temporal_folder}\\subset_{i + 1}.csv', index=False)
            print(f"Already processed {i+1}/{count} files")

    def load_taipei_city_border(self)->gpd.GeoDataFrame:
        """
        Load Taipei city border data
        """
        tp_border_path = os.path.join(self.download_folder, "臺北市區界圖_20220915", "G97_A_CADIST_P.shp")
        tp_border = gpd.read_file(tp_border_path, Encoding='utf-8')
        tp_border = tp_border.set_crs("EPSG:3826")
        tp_border = tp_border[["TNAME", "geometry"]]
        return tp_border

    def find_intersecting_hexagons_within_TP(self, temporal_folder, tp_border)->gpd.GeoDataFrame:
        """
        Find intersecting hexagons within Taipei City
        Approximate time: 
        - find_intersecting_hexagons: 1min
        """
        if self.export_tp_all == True:
            tp_hex_5m = find_intersecting_hexagons(temporal_folder, tp_border)

        return tp_hex_5m

    def find_intersecting_hexagons_with_road(self, temporal_folder, tp_border)->gpd.GeoDataFrame:
        """
        Find intersecting hexagons with road data. There's a slight difference of data shape in different ways.
        1. within 5m buffer from road:(2204412, 51)
        2. intersects road polygon:(1856454, 6)

        Approximate time:
        - road_data_preprocess: 20s
        - find_intersecting_hexagons: 13min
        """
        def road_data_preprocess(self):
            """
            Preprocess road data. To avoid issues with using 'intersects' in spatial joins and 
            single dissolved polygons, we opted to pre-process and integrate road datasets 
            before batch matching.

            """
            # Over 8m roads
            above_8m_road_path = os.path.join(self.download_folder, "8mroadup", "8mroadup", "Road.shp")
            above_8m_road = gpd.read_file(above_8m_road_path)

            # Under 8m roads
            under_8m_road_path = os.path.join(self.download_folder, "TP_Road_under_8m.gpkg")
            under_8m_road = gpd.read_file(under_8m_road_path)
            
            # 
            # Drop useless columns & rename columns
            keep_col_a = ['RoadWidth', 'RoadName', 'Road_ID', 'geometry']
            rename_col = {"RoadWidth": "width", "RoadName": "road_name", "Road_ID": "road_id"}
            above_8m_road = above_8m_road[keep_col_a]
            above_8m_road.rename(columns=rename_col, inplace=True)

            keep_col_u = ['ROADID', 'ROADNAME', 'WIDTH', 'geometry']
            rename_col = {"WIDTH": "width", "ROADNAME": "road_name", "ROADID": "road_id"}
            under_8m_road = under_8m_road[keep_col_u]
            under_8m_road.rename(columns=rename_col, inplace=True)

            # Concatenate above_8m_road & under_8m_road
            road_all = pd.concat([above_8m_road, under_8m_road], ignore_index=True)

            return road_all
        
        # Find intersecting hexagons with road data
        road_all = road_data_preprocess(self)
        tp_road_hex_5m = find_intersecting_hexagons(temporal_folder, road_all)
        tp_road_hex_5m = tp_road_hex_5m.drop(columns=['index_right'])
        tp_road_hex_5m.drop_duplicates(subset='id', inplace=True)

        # Give hexagon district attribute
        tp_dist_hex_5m = gpd.sjoin(tp_road_hex_5m, tp_border, how='inner', predicate='intersects')
        tp_dist_hex_5m.drop_duplicates(subset=['id'], inplace=True)
        tp_dist_hex_5m.drop(columns=['index_right'], inplace=True)

        # Convert centroid to WKT format
        tp_dist_hex_5m["centroid"] = tp_dist_hex_5m["centroid"].apply(wkt.loads)
        tp_dist_hex_5m["centroid"] = tp_dist_hex_5m["centroid"].apply(lambda x: x.wkt)
        col_order = ['id', 'TNAME', 'road_id', 'road_name', 'width', 'centroid','geometry']
        tp_dist_hex_5m = tp_dist_hex_5m[col_order]

        return tp_dist_hex_5m
    
    def filter_hexagons_by_town(self)->None:
        if self.export_tp_road_selected:
            tp_dist_hex_5m = tp_dist_hex_5m[tp_dist_hex_5m["TNAME"].isin(self.select_town)]
            tp_dist_hex_5m = tp_dist_hex_5m.reset_index(drop=True)
            print(f"Selected town is:{self.select_town}\ndata shape is:{tp_dist_hex_5m.shape}")
            # Export to csv
            export_path_csv = os.path.join(self.export_folder, f"tp_road_hex_town_selected.csv")
            tp_dist_hex_5m.to_csv(export_path_csv, index=False)
        else:
            print("No town selected, already export all data in previous step.")
    
    def remove_temporal_folder(self)->None:
        """
        Remove temporal folder
        """
        shutil.rmtree(self.temporal_folder)


    def run(self):
        # Check if the temporal folder exists
        self.check_temporal_folder()

        # Check if the original hexagon file exists
        self.check_original_hex_file()

        # Check memory and storage
        if not self.check_memory_and_storage():
            return
        
        # Load hexagon data
        print("===Loading hexagon data===")
        print("===It takes about 15 minutes===")
        hex_5m = self.load_hex()

        # Split hexagon data
        print("===Splitting hexagon data===")
        print("===It takes about 34 minutes===")
        self.split_hex(hex_5m)

        # Load Taipei city border data
        tp_border = self.load_taipei_city_border()

        # Find intersecting hexagons within Taipei City
        if self.export_tp_all == True:
            print("===Processing all hexagons covering Taipei City.===")
            print("===It takes about 40 minutes.===")
            # Find intersecting hexagons within Taipei City
            tp_hex_5m = self.find_intersecting_hexagons_within_TP(self.temporal_folder, tp_border)
            export_path = os.path.join(self.export_folder, "tp_hex_coverwholecity.csv")
            tp_hex_5m.to_csv(export_path, index=False)

        # Find intersecting hexagons with road data
        print("===Processing all hexagons intersected with roads inside Taipei City.===")
        print("===It takes about 14 minutes===")
        tp_dist_hex_5m = self.find_intersecting_hexagons_with_road(self.temporal_folder, tp_border)
        if self.export_tp_road_all == True:
            export_path = os.path.join(self.export_folder, "tp_road_hex.csv")
            tp_dist_hex_5m.to_csv(export_path, index=False)
        
        # Filter hexagons by town
        self.filter_hexagons_by_town()

        # Remove temporal folder
        print("===Removing temporal folder===")
        self.remove_temporal_folder()

        print("Preprocessing completed.")

if __name__ == "__main__":
    preprocess = PreprocessHex()

    START_TIME = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f'Data Preprocessing Start at {START_TIME}!\n')
    print("===This process will take about 1 hour, please be patient.===")

    preprocess.run()

    END_TIME = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f'\nData Preprocessing Finished at {END_TIME}!')