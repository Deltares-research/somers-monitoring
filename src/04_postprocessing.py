import configparser
import psycopg2
import pandas as pd
from shapely import wkb, wkt
from somers import to_linux
import os

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

#%% Database connection
print("Starting initial connection to PostgreSQL")      

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

monitoringsjaar = 2016
#%% make table

sql_create_table = f"""CREATE TABLE IF NOT EXISTS d_post.d_{monitoringsjaar}(
    parcel_id TEXT UNIQUE NOT NULL,
    minimum FLOAT NOT NULL,
    median FLOAT NOT NULL,
    maximum FLOAT NOT NULL,
    perc_10 FLOAT NOT NULL,
    perc_90 FLOAT NOT NULL,
    measure TEXT,
    geometry GEOMETRY);

    ALTER TABLE d_post.d_{monitoringsjaar} OWNER TO somers"""
 
cursor.execute(sql_create_table)
connection.commit()

#%% load input

sql_select_query = f"""SELECT parcel_id, minimum, median, maximum, perc_10, 
    perc_90, measure, geometry, ST_AsText(geometry) FROM c_output.c_{monitoringsjaar}"""

cursor.execute(sql_select_query)
c_monitoringsjaar = cursor.fetchall()

columns_c = {'parcel_id':str,
            'minimum' : float,
            'median' : float,
            'maximum' : float,
            'perc_10' : float,
            'perc_90' : float,
            'measure': str,
            'geometry': object,
            'geometry_wkt':object,
            }

c_df = pd.DataFrame(c_monitoringsjaar, columns = columns_c.keys()).astype(columns_c).set_index('parcel_id')



sql_select_query = f""" SELECT parcel_id, minimum, median, maximum, perc_10, 
    perc_90, measure, geometry, ST_AsText(geometry) FROM c_output.c_{monitoringsjaar}_small_archetypes"""

cursor.execute(sql_select_query)
c_monitoringsjaar_small = cursor.fetchall()

c_small_df = pd.DataFrame(c_monitoringsjaar_small, columns = columns_c.keys()).astype(columns_c)




sql_select_query = f"""SELECT perceel_id, dekking_ve,dekking_mo FROM a_input.a_{monitoringsjaar}"""

cursor.execute(sql_select_query)
a_monitoringsjaar = cursor.fetchall()

columns_a = {'perceel_id': str,
          'dekking_ve': float,
          'dekking_mo': float,
          }

a_df = pd.DataFrame(a_monitoringsjaar, columns = columns_a.keys()).astype(columns_a).set_index('perceel_id')

#%%Check which parcels have already been calculated
sql_completed_parcels = f"""SELECT DISTINCT parcel_id 
       FROM d_post.d_{monitoringsjaar}"""
       
cursor.execute(sql_completed_parcels)
completed_parcels = [str(i[0]) for i in cursor.fetchall()]
    
print("Initial PostgreSQL connection is closed") 

#%% post processing

# c_df = c_df[~c_df.index.isin(completed_parcels)]
c_df = c_df.merge(a_df[['dekking_ve']], left_index=True, right_index=True, how='left')
c_df = c_df.merge(a_df[['dekking_mo']], left_index=True, right_index=True, how='left')
c_df['dekkingsgraad_max'] = c_df[['dekking_ve', 'dekking_mo']].max(axis=1) / 100
c_df['dekkingsgraad_min'] = c_df[['dekking_ve', 'dekking_mo']].min(axis=1) / 100

#identical_rows = c_df[c_df.duplicated(keep=False)]
# c_df = c_df.drop_duplicates()

double_lst = []
single_lst = []

for parcel_id, row in c_df.iterrows():
    
    print(parcel_id)
    parcel_id = parcel_id
    geo_trans = wkt.loads(row.geometry_wkt) 
    area = geo_trans.area/10000 #ha

    if parcel_id in c_small_df.parcel_id.values:
        double_lst.append(parcel_id)
        emission_min = row['minimum'] * area * row.dekkingsgraad_max + c_small_df.loc[c_small_df['parcel_id'] == parcel_id, 'minimum'].values[0] * area * row.dekkingsgraad_min #kg co2 / year
        emission_mean = row['median']  * area *  row.dekkingsgraad_max + c_small_df.loc[c_small_df['parcel_id'] == parcel_id, 'median'].values[0] * area * row.dekkingsgraad_min 
        emission_max = row['maximum'] * area *  row.dekkingsgraad_max + c_small_df.loc[c_small_df['parcel_id'] == parcel_id, 'maximum'].values[0] * area * row.dekkingsgraad_min 
        emission_10_perc = row['perc_10'] * area *  row.dekkingsgraad_max + c_small_df.loc[c_small_df['parcel_id'] == parcel_id, 'perc_10'].values[0] * area * row.dekkingsgraad_min 
        emission_90_perc = row['perc_90'] * area *  row.dekkingsgraad_max + c_small_df.loc[c_small_df['parcel_id'] == parcel_id, 'perc_90'].values[0] * area * row.dekkingsgraad_min 
    
    else: 
        single_lst.append(parcel_id)
        emission_min = row['minimum'] * area * row.dekkingsgraad_max
        emission_mean = row['median'] * area *  row.dekkingsgraad_max
        emission_max = row['maximum'] * area *  row.dekkingsgraad_max
        emission_10_perc = row['perc_10'] * area *  row.dekkingsgraad_max
        emission_90_perc = row['perc_90'] * area *  row.dekkingsgraad_max 
    
    postgres_insert_query = f"""
        INSERT INTO d_post.d_{monitoringsjaar}
        (parcel_id, minimum, median, maximum, perc_10, perc_90, measure, geometry) 
        VALUES (%(parcel_id)s, %(minimum)s, %(median)s, %(maximum)s, %(perc_10)s, %(perc_90)s, %(measure)s, %(geometry)s) 
        ON CONFLICT (parcel_id) 
        DO UPDATE SET 
            minimum = EXCLUDED.minimum,
            median = EXCLUDED.median,
            maximum = EXCLUDED.maximum,
            perc_10 = EXCLUDED.perc_10,
            perc_90 = EXCLUDED.perc_90,
            measure = EXCLUDED.measure,
            geometry = EXCLUDED.geometry
    """
    record_to_insert = {'parcel_id':str(parcel_id), 'minimum':float(emission_min), 
                        'median':float(emission_mean), 'maximum':float(emission_max),'perc_10':float(emission_10_perc),'perc_90':float(emission_90_perc), 
                        'measure':row.measure,'geometry':row.geometry}
    cursor.execute(postgres_insert_query, record_to_insert)
    
    # output is in kg co2/perceel/jaar
    connection.commit() 
    count = cursor.rowcount
    print(count, "Record inserted successfully into mobile table")

cursor.close()
connection.close()

