import pandas as pd
import warnings
import configparser
import geopandas as gpd
import psycopg2
import sys
import os
import numpy as np
os.environ['USE_PYGEOS'] = '0'
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
warnings.filterwarnings("ignore", category=DeprecationWarning)

from somers import to_linux
sys.path.append(to_linux(r'P:/11207812-somers-uitvoering/monitoring_2023/src/'))
from output_functions import (adjust_freeboard,
                                run_somers_ref,
                                run_somers_ssi,
                                run_somers_pssi,
                                )

#%% Database connection
print("Starting initial connection to PostgreSQL")      

#Read the variables from the ssh command-line arguments
variables = sys.argv[1:]
# variables = [2023, 0, 2000] ## DEZE MOET JE COMMENTEN NA HET PROBEREN
#print(sys.argv[1:])
monitorings_jaar_idx = int(variables[0])
start_idx = int(variables[1])
end_idx = int(variables[2])

# Establish database connection & make table
user_credentials = 'ln'
config_file_path = to_linux(f'P:/11207812-somers-uitvoering/monitoring_2023/configfile_somers_db_{user_credentials}.txt')

#%%

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
#%%

monitoringsjaar = monitorings_jaar_idx
#%%

# sql_create_table = f"""CREATE TABLE IF NOT EXISTS c_output.c_{monitoringsjaar}_gr(
#     parcel_id TEXT UNIQUE NOT NULL,
#     minimum FLOAT NOT NULL,
#     median FLOAT NOT NULL,
#     maximum FLOAT NOT NULL,
#     perc_10 FLOAT NOT NULL,
#     perc_90 FLOAT NOT NULL,
#     min_summer_gwl FLOAT NOT NULL,
#     min_summer3_gwl FLOAT NOT NULL,
#     min_winter_gwl FLOAT NOT NULL,
#     min_year_gwl FLOAT NOT NULL,
#     min_rlg_gwl FLOAT NOT NULL,
#     min_rhg_gwl FLOAT NOT NULL,
#     mean_summer_gwl FLOAT NOT NULL,
#     mean_summer3_gwl FLOAT NOT NULL,
#     mean_winter_gwl FLOAT NOT NULL,
#     mean_year_gwl FLOAT NOT NULL,
#     mean_rlg_gwl FLOAT NOT NULL,
#     mean_rhg_gwl FLOAT NOT NULL,
#     max_summer_gwl FLOAT NOT NULL,
#     max_summer3_gwl FLOAT NOT NULL,
#     max_winter_gwl FLOAT NOT NULL,
#     max_year_gwl FLOAT NOT NULL,
#     max_rlg_gwl FLOAT NOT NULL,
#     max_rhg_gwl FLOAT NOT NULL,
#     succesfull_runs INT NOT NULL,
#     measure TEXT);

#     ALTER TABLE c_output.c_{monitoringsjaar}_gr OWNER TO somers"""

# cursor.execute(sql_create_table)
# connection.commit()

b_monitoringsjaar = gpd.read_file(to_linux(r'P:/11207812-somers-uitvoering/monitoring_2023/src/b_2016_gr.shp'))
b_monitoringsjaar = b_monitoringsjaar.set_index('PARCEL_ID')

#%%Check which parcels have already been calculated

sql_completed_parcels = f"""SELECT DISTINCT parcel_id 
       FROM c_output.c_{monitoringsjaar}_gr"""
cursor.execute(sql_completed_parcels)
completed_parcels = [str(i[0]) for i in cursor.fetchall()]

cursor.close()
connection.close()

ssi_depth_mean = b_monitoringsjaar.SSI_DEPTH.mean()
    
print("Initial PostgreSQL connection is closed") 
print("Starting model run connection to PostgreSQL")

parcels = b_monitoringsjaar[~b_monitoringsjaar.index.isin(completed_parcels)]

skip_lst = ['186868_2']

parcels = parcels[~parcels.index.isin(skip_lst)]
parcels = parcels.iloc[start_idx:end_idx + 1]
#%%
print(f'running n={len(parcels)} parcels')


runs = pd.DataFrame()    
for parcel_id, parcel in parcels.iterrows():   
    print(parcel_id)
    
    #make connection
    connection = psycopg2.connect(host=host, database=database, user=user, password=password, port=port)
    cursor = connection.cursor()
    
    ## load input
    x = parcel.X
    y = parcel.Y
    archetype = parcel.ARCHETYPE
    parcel_width = parcel.PARCEL_WID
    surface_level = parcel.SURFACE_LE #AHN3
    summer_stage = adjust_freeboard(surface_level, parcel.SUMMER_STA)
    winter_stage = adjust_freeboard(surface_level, parcel.WINTER_STA)
    ssi_distance = parcel.SSI_DISTAN
    ssi_depth = parcel.SSI_DEPTH
    pssi_distance = parcel.PSSI_DISTA
    pssi_depth = parcel.PSSI_DEPTH
    measure = parcel.MEASURE
    geometry = parcel.geometry
        
    
    if parcel.MEASURE =='SSI':
        print('Run SSI model sim')
        # continue
    
        if np.isnan(parcel.ssi_depth):
            ssi_depth = ssi_depth_mean
        else:
            ssi_depth = ssi_depth
        
        output_min, output_median, output_max, output_10_perc, output_90_perc,min_summer_gwl,min_summer3_gwl,min_winter_gwl,min_year_gwl,min_rlg_gwl,min_rhg_gwl,mean_summer_gwl,mean_summer3_gwl,mean_winter_gwl,mean_year_gwl,mean_rlg_gwl,mean_rhg_gwl,max_summer_gwl,max_summer3_gwl,max_winter_gwl,max_year_gwl,max_rlg_gwl,max_rhg_gwl, nb_runs = run_somers_ssi(
            parcel_id=parcel_id,
            x=x,
            y=y,
            archetype=archetype,
            parcel_width=parcel_width,
            surface_level=surface_level,
            summer_stage=summer_stage,
            winter_stage=winter_stage,
            ssi_distance=ssi_distance,
            ssi_depth=ssi_depth)

    elif parcel.MEASURE =='PSSI':
        print('Run PSSI model sim')
        # continue      
        output_min, output_median, output_max, output_10_perc, output_90_perc,min_summer_gwl,min_summer3_gwl,min_winter_gwl,min_year_gwl,min_rlg_gwl,min_rhg_gwl,mean_summer_gwl,mean_summer3_gwl,mean_winter_gwl,mean_year_gwl,mean_rlg_gwl,mean_rhg_gwl,max_summer_gwl,max_summer3_gwl,max_winter_gwl,max_year_gwl,max_rlg_gwl,max_rhg_gwl, nb_runs = run_somers_pssi(
            parcel_id=parcel_id,
            x=x,
            y=y,
            archetype=archetype,
            parcel_width=parcel_width,
            surface_level=surface_level,
            summer_stage=summer_stage,
            winter_stage=winter_stage,
            pssi_distance=pssi_distance,
            pssi_depth=pssi_depth)
        
    elif parcel.MEASURE == 'ref':
        print('Run ref model sim')
        #continue
        output_min, output_median, output_max, output_10_perc, output_90_perc,min_summer_gwl,min_summer3_gwl,min_winter_gwl,min_year_gwl,min_rlg_gwl,min_rhg_gwl,mean_summer_gwl,mean_summer3_gwl,mean_winter_gwl,mean_year_gwl,mean_rlg_gwl,mean_rhg_gwl,max_summer_gwl,max_summer3_gwl,max_winter_gwl,max_year_gwl,max_rlg_gwl,max_rhg_gwl, nb_runs = run_somers_ref(
            parcel_id=parcel_id,
            x=x,
            y=y,
            archetype=archetype,
            parcel_width=parcel_width,
            surface_level=surface_level,
            summer_stage=summer_stage,
            winter_stage=winter_stage)
    
    if output_min is None:
        continue
    

    postgres_insert_query = f"""INSERT INTO c_output.c_{monitoringsjaar}_gr (parcel_id, minimum, median, maximum, perc_10, perc_90, 
        min_summer_gwl,min_summer3_gwl,min_winter_gwl,min_year_gwl,min_rlg_gwl,min_rhg_gwl,
        mean_summer_gwl,mean_summer3_gwl,mean_winter_gwl,mean_year_gwl,mean_rlg_gwl,mean_rhg_gwl,
        max_summer_gwl,max_summer3_gwl,max_winter_gwl,max_year_gwl,max_rlg_gwl,max_rhg_gwl,
        succesfull_runs, measure)
        VALUES (%(parcel_id)s, %(minimum)s, %(median)s, %(maximum)s, %(perc_10)s, %(perc_90)s, 
        %(min_summer_gwl)s, %(min_summer3_gwl)s, %(min_winter_gwl)s, %(min_year_gwl)s, %(min_rlg_gwl)s, %(min_rhg_gwl)s, 
        %(mean_summer_gwl)s, %(mean_summer3_gwl)s, %(mean_winter_gwl)s, %(mean_year_gwl)s, %(mean_rlg_gwl)s, %(mean_rhg_gwl)s, 
        %(max_summer_gwl)s, %(max_summer3_gwl)s, %(max_winter_gwl)s, %(max_year_gwl)s, %(max_rlg_gwl)s, %(max_rhg_gwl)s,
        %(succesfull_runs)s, %(measure)s)
        ON CONFLICT (parcel_id) DO UPDATE 
        SET minimum = EXCLUDED.minimum, median = EXCLUDED.median, maximum = EXCLUDED.maximum, perc_10 = EXCLUDED.perc_10, perc_90 = EXCLUDED.perc_90, 
        min_summer_gwl = Excluded.min_summer_gwl, min_summer3_gwl = Excluded.min_summer3_gwl, min_winter_gwl = Excluded.min_winter_gwl, min_year_gwl = Excluded.min_year_gwl, min_rlg_gwl = Excluded.min_rlg_gwl, min_rhg_gwl = Excluded.min_rhg_gwl,
        mean_summer_gwl = Excluded.mean_summer_gwl, mean_summer3_gwl = Excluded.mean_summer3_gwl, mean_winter_gwl = Excluded.mean_winter_gwl, mean_year_gwl = Excluded.mean_year_gwl, mean_rlg_gwl = Excluded.mean_rlg_gwl, mean_rhg_gwl = Excluded.mean_rhg_gwl,
        max_summer_gwl = Excluded.max_summer_gwl, max_summer3_gwl = Excluded.max_summer3_gwl, max_winter_gwl = Excluded.max_winter_gwl, max_year_gwl = Excluded.max_year_gwl, max_rlg_gwl = Excluded.max_rlg_gwl, max_rhg_gwl = Excluded.max_rhg_gwl,
        succesfull_runs = EXCLUDED.succesfull_runs, measure = EXCLUDED.measure"""
    record_to_insert = {'parcel_id': str(parcel_id), 'minimum': float(output_min), 'median': float(output_median), 'maximum': float(output_max),'perc_10': float(output_10_perc), 'perc_90': float(output_90_perc),
                        'min_summer_gwl': float(min_summer_gwl), 'min_summer3_gwl': float(min_summer3_gwl), 'min_winter_gwl': float(min_winter_gwl), 'min_year_gwl': float(min_year_gwl), 'min_rlg_gwl': float(min_rlg_gwl), 'min_rhg_gwl': float(min_rhg_gwl),
                        'mean_summer_gwl': float(mean_summer_gwl), 'mean_summer3_gwl': float(mean_summer3_gwl), 'mean_winter_gwl': float(mean_winter_gwl), 'mean_year_gwl': float(mean_year_gwl), 'mean_rlg_gwl': float(mean_rlg_gwl), 'mean_rhg_gwl': float(mean_rhg_gwl),
                        'max_summer_gwl': float(max_summer_gwl), 'max_summer3_gwl': float(max_summer3_gwl), 'max_winter_gwl': float(max_winter_gwl), 'max_year_gwl': float(max_year_gwl), 'max_rlg_gwl': float(max_rlg_gwl), 'max_rhg_gwl': float(max_rhg_gwl),
                        'succesfull_runs': int(nb_runs), 'measure': str(measure)}
    cursor.execute(postgres_insert_query, record_to_insert)

    # output is in kg co2/ha/jaar
    connection.commit()
    count = cursor.rowcount
    print(count, "Record inserted successfully into mobile table")

    cursor.close()
    connection.close()
    print("PostgreSQL connection for single run is closed") 

print("PostgreSQL connection for all runs is closed") 

#%%

# sql_select_query = f"""SELECT parcel_id, minimum, median, maximum, perc_10, 
#     perc_90, measure FROM c_output.c_{monitoringsjaar}_gr"""

# cursor.execute(sql_select_query)
# c_monitoringsjaar = cursor.fetchall()

# columns_c = {'parcel_id':str,
#             'minimum' : float,
#             'median' : float,
#             'maximum' : float,
#             'perc_10' : float,
#             'perc_90' : float,
#             'measure': str,
#             }

# c_df = pd.DataFrame(c_monitoringsjaar, columns = columns_c.keys()).astype(columns_c).set_index('parcel_id')

# b_monitoringsjaar = gpd.read_file(to_linux(r'P:/11207812-somers-uitvoering/monitoring_2023/src/b_2016_fr.shp'))
# b_monitoringsjaar = b_monitoringsjaar.set_index('PARCEL_ID')

# # b_monitoringsjaar = gpd.read_file(to_linux(r'P:/11207812-somers-uitvoering/monitoring_2021/src/b_2016_fr.shp'))
# # b_monitoringsjaar = b_monitoringsjaar.set_index('PARCEL_ID')


# c_df['geometry'] = b_monitoringsjaar['geometry']

# a_df = gpd.read_file(to_linux(r'P:/11207812-somers-uitvoering/monitoring_2021/database/2-intermin/a_input/input_parcels_2016.shp'))
# a_df = a_df.set_index('Perceel_ID')

# # a_df = gpd.read_file(to_linux(r'P:/11207812-somers-uitvoering/monitoring_2023/shp_files/a_input/a_2016.shp'))
# # a_df = a_df.set_index('PERCEEL_ID')

# c_df['dekking_veen'] = a_df['Dekkingsgr'] 
# c_df['dekking_moerig'] = a_df['Dekkings_1']


# # c_df['dekking_veen'] = a_df['DEKKING_VE'] 
# # c_df['dekking_moerig'] = a_df['DEKKING_MO']

# c_gdf = gpd.GeoDataFrame(c_df, geometry='geometry')

# c_gdf['org_opp'] = (c_gdf.geometry.area * ((c_gdf['dekking_veen'] + c_gdf['dekking_moerig']) / 100)) / 10000 #ha
# c_gdf['uitstoot_tot'] = (c_gdf['median'] * c_gdf['org_opp']) / 1000 #ton CO2

# tot_uitstoot = c_gdf['uitstoot_tot'].sum() / 1000 #k ton CO2

# gem_uitstoot = (np.sum(c_gdf['median'] * c_gdf['org_opp']) / np.sum(c_gdf['org_opp'])) / 1000 #ton CO2 / ha

# org_area = np.sum(c_gdf['org_opp']) #ha

# #%%

# c_gdf['org_dekking'] = ((c_gdf['dekking_veen'] + c_gdf['dekking_moerig']) / 100)
# check = c_gdf['org_dekking'].mean()
