import configparser
import psycopg2
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from shapely import wkb, wkt
from shapely.wkt import loads
import geopandas as gpd
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

monitoringsjaar = 2023
#%%

sql_select_query = f"""SELECT parcel_id, minimum, median, maximum, perc_10, 
    perc_90, measure, geometry, ST_AsText(geometry) FROM d_post.d_{monitoringsjaar}_ahn3"""

cursor.execute(sql_select_query)
d_monitoringsjaar = cursor.fetchall()

columns_d = {'parcel_id':str,
            'minimum' : float,
            'median' : float,
            'maximum' : float,
            'perc_10' : float,
            'perc_90' : float,
            'measure': str,
            'geometry': object,
            'geometry_wkt':object,
            }

d_df = pd.DataFrame(d_monitoringsjaar, columns = columns_d.keys()).astype(columns_d).set_index('parcel_id')
d_df['geometry'] = d_df['geometry_wkt'].apply(loads)
d_df = gpd.GeoDataFrame(d_df, geometry='geometry', crs="EPSG:28992")  # Assuming WGS84 (latitude/longitude)
d_df.drop(columns=['geometry_wkt'], inplace=True)



sql_select_query = f"""SELECT perceel_id, dekking_ve,dekking_mo FROM a_input.a_{monitoringsjaar}"""

cursor.execute(sql_select_query)
a_monitoringsjaar = cursor.fetchall()

columns_a = {'perceel_id': str,
          'dekking_ve': float,
          'dekking_mo': float,
          }

a_df = pd.DataFrame(a_monitoringsjaar, columns = columns_a.keys()).astype(columns_a).set_index('perceel_id')


#%%
d_df['dekking_ve'] = d_df.index.map(a_df['dekking_ve'])
d_df['dekking_mo'] = d_df.index.map(a_df['dekking_mo'])
  
d_df['dekking_org_totaal'] = (d_df['dekking_ve'] + d_df['dekking_mo']) / 100
d_df['median_per_ha'] = d_df['median'] / ((d_df.geometry.area / 10000) * d_df['dekking_org_totaal'])

d_emission = d_df['median_per_ha'] * (d_df.geometry.area / 10000) * ((d_df['dekking_ve'] / 100) + (d_df['dekking_mo'] / 100))
d_peat_emission = d_df['median_per_ha'] * (d_df.geometry.area / 10000) * (d_df['dekking_ve'] / 100)
d_peaty_emission = d_df['median_per_ha'] * (d_df.geometry.area / 10000) * (d_df['dekking_mo'] / 100)

gewogen_mean =  np.sum(d_df['median_per_ha'] * (d_df.geometry.area / 10000) * d_df['dekking_org_totaal']) / np.sum((d_df.geometry.area / 10000) * d_df['dekking_org_totaal'])
gewogen_mean_peat =  np.sum(d_peat_emission) / np.sum((d_df.geometry.area / 10000) * (d_df['dekking_ve'] / 100))
gewogen_mean_peaty =  np.sum(d_peaty_emission) / np.sum((d_df.geometry.area / 10000) * (d_df['dekking_mo'] / 100))

peat_area = np.sum((d_df.geometry.area / 10000) * (d_df['dekking_ve'] / 100))
peaty_area = np.sum((d_df.geometry.area / 10000) * (d_df['dekking_mo'] / 100))
total_area = np.sum((d_df.geometry.area / 10000) * d_df['dekking_org_totaal'])

#%%
# print values
print(f'Total emissisions in {monitoringsjaar}:', f"{d_df['median'].sum() / (1000000 * 1000) :.2f}", 'Mton CO2')
print(f'Total emissisions in {monitoringsjaar}:', f"{d_emission.sum() / (1000000 * 1000) :.2f}", 'Mton CO2')
print('Of which', f"{d_peat_emission.sum() / (1000000 * 1000):.2f}", 'Mton CO2 is from peat soils')
print('and', f"{d_peaty_emission.sum() / (1000000 * 1000):.2f}", 'Mton CO2 is from peaty soils')
print(f'Min emissisions in {monitoringsjaar}:', f"{d_df['minimum'].sum() / (1000000 * 1000):.2f}", 'Mton CO2')
print(f'Max emissisions in {monitoringsjaar}:', f"{d_df['maximum'].sum() / (1000000 * 1000):.2f}", 'Mton CO2')
print(f'10th perc emissisions in {monitoringsjaar}:', f"{d_df['perc_10'].sum() / (1000000 * 1000):.2f}", 'Mton CO2')
print(f'90th perc emissisions in {monitoringsjaar}:', f"{d_df['perc_90'].sum() / (1000000 * 1000):.2f}", 'Mton CO2')

print(f'Emissions per hectare in {monitoringsjaar}:', f"{gewogen_mean / 1000:.2f}", 'ton CO2/ha')
print('where specifically', f"{gewogen_mean_peat / 1000:.2f}",'ton CO2/ha for peat soils')
print('and', f"{gewogen_mean_peaty/ 1000:.2f}", 'ton CO2/ha for peaty soils')

print('Organic area is', total_area, 'ha')
print('of which', peat_area, 'ha falls on peat soils')
print('and', peaty_area, 'ha falls on peaty soils')
#%% general per provincie

provincies = gpd.read_file(r'P:/11207812-somers-ontwikkeling/base_data/GIS_layers/Bestuurlijke_grenzen/B1_Provinciale_indeling_van_NederlandPolygon.shp')

provincie_ls = []

for parcel_id, row in d_df.iterrows():
    #print (row)
    geom_row = row.geometry.centroid
    province = provincies[provincies.contains(geom_row)]['PROVC_NM']
    if province.empty:
        province = np.nan
    else:
        province = province.iloc[0]
    provincie_ls.append(province)

d_df['provincie'] = provincie_ls

#%%

d_zuid_holland = d_df[d_df['provincie'] == 'Zuid-Holland']
d_utrecht = d_df[d_df['provincie'] == 'Utrecht']
d_noord_holland = d_df[d_df['provincie'] == 'Noord-Holland']
d_overijssel = d_df[d_df['provincie'] == 'Overijssel']
d_drenthe = d_df[d_df['provincie'] == 'Drenthe']
d_groningen = d_df[d_df['provincie'] == 'Groningen']
d_friesland = d_df[d_df['provincie'] == 'Friesland']
d_gelderland = d_df[d_df['provincie'] == 'Gelderland']
d_zeeland = d_df[d_df['provincie'] == 'Zeeland']
d_flevoland = d_df[d_df['provincie'] == 'Flevoland']
d_brabant = d_df[d_df['provincie'] == 'Noord-Brabant']


gewogen_mean_groningen =  np.sum(d_groningen['median_per_ha'] * d_groningen.geometry.area / 10000 * d_groningen['dekking_org_totaal']) / np.sum(d_groningen.geometry.area / 10000 * d_groningen['dekking_org_totaal'])
gewogen_mean_friesland =  np.sum(d_friesland['median_per_ha'] * d_friesland.geometry.area / 10000 * d_friesland['dekking_org_totaal']) / np.sum(d_friesland.geometry.area / 10000 * d_friesland['dekking_org_totaal'])
gewogen_mean_overijssel =  np.sum(d_overijssel['median_per_ha'] * d_overijssel.geometry.area / 10000 * d_overijssel['dekking_org_totaal']) / np.sum(d_overijssel.geometry.area / 10000 * d_overijssel['dekking_org_totaal'])
gewogen_mean_drenthe =  np.sum(d_drenthe['median_per_ha'] * d_drenthe.geometry.area / 10000 * d_drenthe['dekking_org_totaal']) / np.sum(d_drenthe.geometry.area / 10000 * d_drenthe['dekking_org_totaal'])
gewogen_mean_gelderland =  np.sum(d_gelderland['median_per_ha'] * d_gelderland.geometry.area / 10000 * d_gelderland['dekking_org_totaal']) / np.sum(d_gelderland.geometry.area / 10000 * d_gelderland['dekking_org_totaal'])
gewogen_mean_utrecht =  np.sum(d_utrecht['median_per_ha'] * d_utrecht.geometry.area / 10000 * d_utrecht['dekking_org_totaal']) / np.sum(d_utrecht.geometry.area / 10000 * d_utrecht['dekking_org_totaal'])
gewogen_mean_noord_holland =  np.sum(d_noord_holland['median_per_ha'] * d_noord_holland.geometry.area / 10000 * d_noord_holland['dekking_org_totaal']) / np.sum(d_noord_holland.geometry.area / 10000 * d_noord_holland['dekking_org_totaal'])
gewogen_mean_zuid_holland =  np.sum(d_zuid_holland['median_per_ha'] * d_zuid_holland.geometry.area / 10000 * d_zuid_holland['dekking_org_totaal']) / np.sum(d_zuid_holland.geometry.area / 10000 * d_zuid_holland['dekking_org_totaal'])
gewogen_mean_zeeland =  np.sum(d_zeeland['median_per_ha'] * d_zeeland.geometry.area / 10000 * d_zeeland['dekking_org_totaal']) / np.sum(d_zeeland.geometry.area / 10000 * d_zeeland['dekking_org_totaal'])
gewogen_mean_flevoland =  np.sum(d_flevoland['median_per_ha'] * d_flevoland.geometry.area / 10000 * d_flevoland['dekking_org_totaal']) / np.sum(d_flevoland.geometry.area / 10000 * d_flevoland['dekking_org_totaal'])
gewogen_mean_brabant =  np.sum(d_brabant['median_per_ha'] * d_brabant.geometry.area / 10000 * d_brabant['dekking_org_totaal']) / np.sum(d_brabant.geometry.area / 10000 * d_brabant['dekking_org_totaal'])


#%%
print(f'Emissions per hectare in {monitoringsjaar} for Groningen:', f"{gewogen_mean_groningen / 1000:.2f}", 'ton CO2/ha')
print(f'Emissions per hectare in {monitoringsjaar} for Friesland:', f"{gewogen_mean_friesland / 1000:.2f}", 'ton CO2/ha')
print(f'Emissions per hectare in {monitoringsjaar} for Overijssel:', f"{gewogen_mean_overijssel / 1000:.2f}", 'ton CO2/ha')
print(f'Emissions per hectare in {monitoringsjaar} for Drenthe:', f"{gewogen_mean_drenthe / 1000:.2f}", 'ton CO2/ha')
print(f'Emissions per hectare in {monitoringsjaar} for Gelderland:', f"{gewogen_mean_gelderland / 1000:.2f}", 'ton CO2/ha')
print(f'Emissions per hectare in {monitoringsjaar} for Utrecht:', f"{gewogen_mean_utrecht / 1000:.2f}", 'ton CO2/ha')
print(f'Emissions per hectare in {monitoringsjaar} for Noord Holland:', f"{gewogen_mean_noord_holland / 1000:.2f}", 'ton CO2/ha')
print(f'Emissions per hectare in {monitoringsjaar} for Zeeland:', f"{gewogen_mean_zeeland / 1000:.2f}", 'ton CO2/ha')
print(f'Emissions per hectare in {monitoringsjaar} for Zuid Holland:', f"{gewogen_mean_zuid_holland / 1000:.2f}", 'ton CO2/ha')
print(f'Emissions per hectare in {monitoringsjaar} for Flevoland:', f"{gewogen_mean_flevoland / 1000:.2f}", 'ton CO2/ha')
print(f'Emissions per hectare in {monitoringsjaar} for Brabant:', f"{gewogen_mean_brabant / 1000:.2f}", 'ton CO2/ha')


print('Total emissisions in groningen:', f"{d_groningen['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in friesland:', f"{d_friesland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in overijssel:', f"{d_overijssel['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in drenthe:', f"{d_drenthe['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in gelderland:', f"{d_gelderland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in utrecht:', f"{d_utrecht['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in noord_holland:', f"{d_noord_holland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in zeeland:', f"{d_zeeland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in zuid_holland:', f"{d_zuid_holland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in flevoland:', f"{d_flevoland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in brabant:', f"{d_brabant['median'].sum() / (1000):.0f}", 'ton CO2')

