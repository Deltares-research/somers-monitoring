import pandas as pd
import warnings
import configparser
import psycopg2
import sys
import os
import numpy as np
os.environ['USE_PYGEOS'] = '0'
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
warnings.filterwarnings("ignore", category=DeprecationWarning)

import tempfile  # For creating a temporary directory
from pathlib import Path

import geopandas as gpd  # For reading the input parcels
import pandas as pd  # For specifying dates of the modelling period
import shapely
from shapely.geometry import Point
import xarray as xr

import somers
from output_functions_refactor import (calculate_stats)

#%% Read SH variables
print("Starting initial connection to PostgreSQL")      

#Read the variables from the command-line arguments
variables = sys.argv[1:]
variables = [2022, 0, 2000] ## DEZE MOET JE COMMENTEN NA HET PROBEREN
#print(sys.argv[1:])
monitorings_jaar_idx = int(variables[0])
start_idx = int(variables[1])
end_idx = int(variables[2])

#%% Establish database connection
user_credentials = 'ln'
config_file_path = Path(f'P:/11207812-somers-uitvoering/monitoring_2023/configfile_somers_db_{user_credentials}.txt')

# Read the configuration file
cf = configparser.ConfigParser()
cf.read(config_file_path)

# Get the database connection parameters from the configuration file
host = cf.get("PostGIS", 'host')
database = cf.get("PostGIS", 'db')
user = cf.get("PostGIS", 'user')
password = cf.get("PostGIS", 'pass')
port = cf.get("PostGIS", 'port')

# Establish a connection to the PostgreSQL database
connection = psycopg2.connect(host=host, database=database, user=user, password=password, port=port)
cursor = connection.cursor()
#%% Select monitoring year from SH variable

monitoringsjaar = monitorings_jaar_idx
#%% Create new database table and make it accessible to everyone

sql_create_table = f"""CREATE TABLE IF NOT EXISTS c_output_LULUCF.c_{monitoringsjaar}(
    parcel_id TEXT UNIQUE NOT NULL,
    minimum FLOAT NOT NULL,
    median FLOAT NOT NULL,
    maximum FLOAT NOT NULL,
    perc_10 FLOAT NOT NULL,
    perc_90 FLOAT NOT NULL,
    min_summer_gwl FLOAT NOT NULL,
    min_summer3_gwl FLOAT NOT NULL,
    min_winter_gwl FLOAT NOT NULL,
    min_year_gwl FLOAT NOT NULL,
    min_rlg_gwl FLOAT NOT NULL,
    min_rhg_gwl FLOAT NOT NULL,
    mean_summer_gwl FLOAT NOT NULL,
    mean_summer3_gwl FLOAT NOT NULL,
    mean_winter_gwl FLOAT NOT NULL,
    mean_year_gwl FLOAT NOT NULL,
    mean_rlg_gwl FLOAT NOT NULL,
    mean_rhg_gwl FLOAT NOT NULL,
    max_summer_gwl FLOAT NOT NULL,
    max_summer3_gwl FLOAT NOT NULL,
    max_winter_gwl FLOAT NOT NULL,
    max_year_gwl FLOAT NOT NULL,
    max_rlg_gwl FLOAT NOT NULL,
    max_rhg_gwl FLOAT NOT NULL,
    measure TEXT,
    geometry GEOMETRY);

    ALTER TABLE c_output_LULUCF.c_{monitoringsjaar} OWNER TO somers"""

cursor.execute(sql_create_table)
connection.commit()

# Download B scheme
sql_select_query = f"""SELECT parcel_id, x, y, archetype, surface_level_ahn3, surface_level_ahn5, parcel_width, 
    summer_stage, winter_stage, ssi_distance, ssi_depth, pssi_distance, pssi_depth, measure, geometry, ST_AsText(geometry) 
    FROM b_selection.b_{monitoringsjaar} ORDER BY parcel_id"""

cursor.execute(sql_select_query)
record = cursor.fetchall()

connection.commit()

columns = {'parcel_id': str,
          'x': float,
          'y': float,
          'archetype': str,
          'surface_level_ahn3': float,
          'surface_level_ahn5': float,
          'parcel_width': float,
          'summer_stage': float,
          'winter_stage': float,
          'ssi_distance': float,
          'ssi_depth': float,
          'pssi_distance': float,
          'pssi_depth': float,
          'measure': str,
          'geometry': object,
          'geometry_wkt': object,
          }

b_monitoringsjaar = pd.DataFrame(record, columns = columns.keys()).astype(columns)
b_monitoringsjaar['geometry_wkt'] = b_monitoringsjaar['geometry_wkt'].apply(shapely.wkt.loads)
b_monitoringsjaar = gpd.GeoDataFrame(b_monitoringsjaar, geometry='geometry_wkt')
b_monitoringsjaar = b_monitoringsjaar.explode(index_parts=False).reset_index(drop=True)


#%% Check which parcels have already been calculated and remove them from the parcels to run

sql_completed_parcels = f"""SELECT DISTINCT parcel_id 
       FROM c_output_LULUCF.c_{monitoringsjaar}"""
cursor.execute(sql_completed_parcels)
completed_parcels = [str(i[0]) for i in cursor.fetchall()]

cursor.close()
connection.close()
print("Initial PostgreSQL connection is closed") 

ssi_depth_mean = b_monitoringsjaar.ssi_depth.mean()
    
#%% Remove parcels that cause MODFLOW to get stuck

# skip_lst = ['RP-L-3077305', 'RP-L-3089835', 'RP-L-3296055', 'RP-L-3417857', 'RP-L-3305432', 'RP-L-3417857', 'RP-L-3417726',
#             'RP-L-3073591', 'RP-L-3353583', 'RP-L-3453667', 'RP-L-3129115', 'RP-H-2297354', 'RP-L-3305970 '] 

parcels = b_monitoringsjaar[~b_monitoringsjaar.index.isin(completed_parcels)]

# parcels = parcels[~parcels.index.isin(skip_lst)]

parcels = parcels.rename(columns={
    "archetype": 'soil_unit',
    'surface_level_ahn3': 'surface_level', #ahn3
})

parcels = parcels.iloc[start_idx:end_idx + 1]

#%% Import files

path_lhm_data = Path(r'P:\11207812-somers-ontwikkeling\somers_v2.1_input\data\LHM')
path_weather_data = Path(r'P:\11207812-somers-ontwikkeling\somers_v2.1_input\data\meteo')
path_soil_map_data = Path(r'P:\11207812-somers-ontwikkeling\somers_v2.1_input\data\ondergrond')
path_ghg_data = Path(r'P:\11207812-somers-ontwikkeling\somers_v2.1_input\data')
path_wfps_data = Path(r'P:\11207812-somers-ontwikkeling\somers_v2.1_input\data\WFPS')


lhm = somers.io.read_lhm_data(
    path_lhm_data / r"lhm4_confining_layer.nc",
    path_lhm_data / r"lhm4.2_flux_l1_2013_2022.nc",
    path_lhm_data / r"lhm4.2_recharge_2010_2023.nc",)

weather = somers.io.read_weather_data(
    path_weather_data / r"weather_stations.geoparquet",
    path_weather_data / r"temperature_2016-2023.txt",
    #path_weather_data / r"weerstations_2012-2016.txt",
    path_weather_data / r"weather_regions.geoparquet",)

soilmap = somers.io.read_bro_soilmap(path_soil_map_data / r"BRO_DownloadBodemkaart.gpkg") #voor nu nog de nieuwe. Dit moet aangepast worden naar de "oude methode"

ghg = somers.io.read_ghg(path_ghg_data / r"GHG-kk2010.nc")

wfps_lookup = somers.io.read_wfps_lookup_tables(
    southwest=path_wfps_data / r"WFPS_weerregio_southwest.nc",
    northeast=path_wfps_data / "WFPS_weerregio_northeast.nc",
    engine="netcdf4",)

modflow_params_ref = somers.io.read_modflow_parameters('P:/11207812-somers-ontwikkeling/somers_v2.1_development/modflow_calibration/4-output/somers_2_1_parameter_selection.csv')
modflow_params_ssi = somers.io.read_modflow_parameters('P:/11207812-somers-ontwikkeling/somers_v2_development/data/3-input/pp2d_phreatic_head/parameters/definitive parameters/parameters_ssi_sv2.csv')
modflow_params_pssi = somers.io.read_modflow_parameters('P:/11207812-somers-ontwikkeling/somers_v2_development/data/3-input/pp2d_phreatic_head/parameters/definitive parameters/parameters_ssi_sv2.csv')

workdir = Path(r'P:/11207812-somers-uitvoering/monitoring_2023/work_dir_refactor')

modflow_exe_local = Path(r"p:\11207812-somers-uitvoering\software\mf6.exe")  
modflow_exe_linux = Path(r"p:\11207812-somers-uitvoering\software\mf6_compiled_rm_v3") 

# start_date=pd.to_datetime(f'01-01-{find_weerjaar(parcel.x, parcel.y)}',format="%d-%m-%Y")
# end_date=pd.to_datetime(f'31-12-{find_weerjaar(parcel.x, parcel.y)}', format="%d-%m-%Y")

start_date=pd.to_datetime('01-01-2022',format="%d-%m-%Y")
end_date=pd.to_datetime('01-12-2022', format="%d-%m-%Y")


    
#%%
print("Starting model run connection to PostgreSQL")
print(f'running n={len(parcels)} parcels')

settings = somers.ModelSettings(workdir, 
    start_date, 
    end_date, 
    modflow_executable=modflow_exe_local)


for parcel_id in parcels['parcel_id']:   
    print(parcel_id)
    parcel = parcels.loc[parcels['parcel_id'] == parcel_id]

    # Make connection
    connection = psycopg2.connect(host=host, database=database, user=user, password=password, port=port)
    cursor = connection.cursor()
    
    measure = parcel['measure'].item()
    geometry = parcel['geometry_wkt'].item()
    
    try: 
        # Select the modules
        if measure == 'SSI':
            print('Run SSI model sim')
        
            if np.isnan(parcel['ssi_depth'].item()):
                ssi_depth = ssi_depth_mean
            else:
                ssi_depth = parcel['ssi_depth'].item()
            
            parcel['pssi_summer_stage'] = parcel['summer_stage'].item()
            parcel['pssi_winter_stage'] = parcel['winter_stage'].item()
            parcel['drain_depth'] = parcel['ssi_depth'].item()
            parcel['drain_distance'] = parcel['ssi_distance'].item()

            
            model_ssi = somers.SomersModel(
                groundwater=somers.pp2d.Modflow(modflow_params_ssi, "flux", measure='ssi'),
                temperature=somers.pp2d.Fft(spinup_days=10),
                soil_moisture=somers.pp2d.DynamicMoisture(),
                aap=somers.aap.Aap("somers", "somers", "somers_v2"),
                settings=settings)

            results = somers.run_somers(parcel, model_ssi, soilmap, lhm, weather, wfps_lookup, ghg)
            
        elif measure == 'PSSI':
            print('Run PSSI model sim')
            
            parcel['pssi_summer_stage'] = parcel['summer_stage'].item()
            parcel['pssi_winter_stage'] = parcel['winter_stage'].item()
            parcel['drain_depth'] = parcel['pssi_depth'].item()
            parcel['drain_distance'] = parcel['pssi_distance'].item()
            
            model_pssi = somers.SomersModel(
                groundwater=somers.pp2d.Modflow(modflow_params_pssi, "flux", measure='pssi'),
                temperature=somers.pp2d.Fft(spinup_days=10),
                soil_moisture=somers.pp2d.DynamicMoisture(),
                aap=somers.aap.Aap("somers", "somers", "somers_v2"),
                settings=settings)

            results = somers.run_somers(parcel, model_pssi, soilmap, lhm, weather, wfps_lookup, ghg)
            
        elif measure == 'ref':
            print('Run ref model sim')

            model_ref = somers.SomersModel(
                groundwater=somers.pp2d.Modflow(modflow_params_ref, "flux"),
                temperature=somers.pp2d.Fft(spinup_days=10),
                soil_moisture=somers.pp2d.DynamicMoisture(),
                aap=somers.aap.Aap("somers", "somers", "somers_v2"),
                settings=settings)

            results = somers.run_somers(parcel, model_ref, soilmap, lhm, weather, wfps_lookup, ghg)

        # Process results (only one model will be run each iteration)
        output = next(results)  #error met 2022 for some reason... 2015 draait wel. Alle input data wordt op schrijf opgeslagen. Hoe kunnen we dit uitzetten, want dit gaat veel opslag kosten
        
        output_median = output['best_estimate'].sum(['time','depth','x']).median('runs').item()
        output_10_perc = output['best_estimate'].sum(['time', 'depth', 'x']).quantile(0.1, 'runs').item()
        output_90_perc = output['best_estimate'].sum(['time', 'depth', 'x']).quantile(0.9, 'runs').item()
        output_min = output['best_estimate'].sum(['time','depth','x']).median('runs').item()
        output_max = output['best_estimate'].sum(['time','depth','x']).median('runs').item()
        
        min_summer_gwl,min_summer3_gwl,min_winter_gwl,min_year_gwl,min_rlg_gwl,min_rhg_gwl = calculate_stats(float(parcel.surface_level), output['phreatic_head'].min('runs'))
        mean_summer_gwl,mean_summer3_gwl,mean_winter_gwl,mean_year_gwl,mean_rlg_gwl,mean_rhg_gwl = calculate_stats(float(parcel.surface_level), output['phreatic_head'].median('runs'))
        max_summer_gwl,max_summer3_gwl,max_winter_gwl,max_year_gwl,max_rlg_gwl,max_rhg_gwl = calculate_stats(float(parcel.surface_level), output['phreatic_head'].max('runs'))

        
        postgres_insert_query = f"""INSERT INTO c_output_LULUCF.c_{monitoringsjaar} (parcel_id, minimum, median, maximum, perc_10, perc_90, 
             min_summer_gwl,min_summer3_gwl,min_winter_gwl,min_year_gwl,min_rlg_gwl,min_rhg_gwl,
             mean_summer_gwl,mean_summer3_gwl,mean_winter_gwl,mean_year_gwl,mean_rlg_gwl,mean_rhg_gwl,
             max_summer_gwl,max_summer3_gwl,max_winter_gwl,max_year_gwl,max_rlg_gwl,max_rhg_gwl,
             measure, geometry)
             VALUES (%(parcel_id)s, %(minimum)s, %(median)s, %(maximum)s, %(perc_10)s, %(perc_90)s, 
             %(min_summer_gwl)s, %(min_summer3_gwl)s, %(min_winter_gwl)s, %(min_year_gwl)s, %(min_rlg_gwl)s, %(min_rhg_gwl)s, 
             %(mean_summer_gwl)s, %(mean_summer3_gwl)s, %(mean_winter_gwl)s, %(mean_year_gwl)s, %(mean_rlg_gwl)s, %(mean_rhg_gwl)s, 
             %(max_summer_gwl)s, %(max_summer3_gwl)s, %(max_winter_gwl)s, %(max_year_gwl)s, %(max_rlg_gwl)s, %(max_rhg_gwl)s,
             %(measure)s, %(geometry)s)
             ON CONFLICT (parcel_id) DO UPDATE 
             SET minimum = EXCLUDED.minimum, median = EXCLUDED.median, maximum = EXCLUDED.maximum, perc_10 = EXCLUDED.perc_10, perc_90 = EXCLUDED.perc_90, 
             min_summer_gwl = Excluded.min_summer_gwl, min_summer3_gwl = Excluded.min_summer3_gwl, min_winter_gwl = Excluded.min_winter_gwl, min_year_gwl = Excluded.min_year_gwl, min_rlg_gwl = Excluded.min_rlg_gwl, min_rhg_gwl = Excluded.min_rhg_gwl,
             mean_summer_gwl = Excluded.mean_summer_gwl, mean_summer3_gwl = Excluded.mean_summer3_gwl, mean_winter_gwl = Excluded.mean_winter_gwl, mean_year_gwl = Excluded.mean_year_gwl, mean_rlg_gwl = Excluded.mean_rlg_gwl, mean_rhg_gwl = Excluded.mean_rhg_gwl,
             max_summer_gwl = Excluded.max_summer_gwl, max_summer3_gwl = Excluded.max_summer3_gwl, max_winter_gwl = Excluded.max_winter_gwl, max_year_gwl = Excluded.max_year_gwl, max_rlg_gwl = Excluded.max_rlg_gwl, max_rhg_gwl = Excluded.max_rhg_gwl,
             measure = EXCLUDED.measure, geometry = EXCLUDED.geometry"""
         
        record_to_insert = {'parcel_id': str(parcel_id), 'minimum': float(output_min), 'median': float(output_median), 'maximum': float(output_max),'perc_10': float(output_10_perc), 'perc_90': float(output_90_perc),
                             'min_summer_gwl': float(min_summer_gwl), 'min_summer3_gwl': float(min_summer3_gwl), 'min_winter_gwl': float(min_winter_gwl), 'min_year_gwl': float(min_year_gwl), 'min_rlg_gwl': float(min_rlg_gwl), 'min_rhg_gwl': float(min_rhg_gwl),
                             'mean_summer_gwl': float(mean_summer_gwl), 'mean_summer3_gwl': float(mean_summer3_gwl), 'mean_winter_gwl': float(mean_winter_gwl), 'mean_year_gwl': float(mean_year_gwl), 'mean_rlg_gwl': float(mean_rlg_gwl), 'mean_rhg_gwl': float(mean_rhg_gwl),
                             'max_summer_gwl': float(max_summer_gwl), 'max_summer3_gwl': float(max_summer3_gwl), 'max_winter_gwl': float(max_winter_gwl), 'max_year_gwl': float(max_year_gwl), 'max_rlg_gwl': float(max_rlg_gwl), 'max_rhg_gwl': float(max_rhg_gwl),
                             'measure': str(measure), 'geometry': geometry.wkt}
        cursor.execute(postgres_insert_query, record_to_insert)
         
         # output is in kg co2/ha/jaar
        connection.commit()
        count = cursor.rowcount
        print(count, "Record inserted successfully into mobile table")
         
        cursor.close()
        connection.close()
        print("PostgreSQL connection for single run is closed") 
         
    except StopIteration:
        print(f"[Warning] No results for parcel {parcel.parcel_id}.")
    except Exception as e:
        print(f"[Error] Parcel {parcel.parcel_id} failed: {e}")


print("PostgreSQL connection for all runs is closed") 
