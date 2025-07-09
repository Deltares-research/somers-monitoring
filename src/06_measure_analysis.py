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

#%% Select monitoringsjaar

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
d_df = gpd.GeoDataFrame(d_df, geometry='geometry', crs="EPSG:28992") 
d_df.drop(columns=['geometry_wkt'], inplace=True)

sql_select_query = f"""SELECT perceel_id, dekking_ve,dekking_mo, peil_gemid, peil_gem_1 FROM a_input.a_{monitoringsjaar}"""

cursor.execute(sql_select_query)
a_monitoringsjaar = cursor.fetchall()

columns_a = {'perceel_id': str,
          'dekking_ve': float,
          'dekking_mo': float,
          'peil_gemid': float, #gem zomerpeil
          'peil_gem_1': float, #gem winterpeil
          }

a_df = pd.DataFrame(a_monitoringsjaar, columns = columns_a.keys()).astype(columns_a).set_index('perceel_id')

#%% ref jaar

sql_select_query_ref = """SELECT parcel_id, minimum, median, maximum, perc_10, 
    perc_90, measure, geometry, ST_AsText(geometry) FROM d_post.d_2016"""

cursor.execute(sql_select_query_ref)
d_ref = cursor.fetchall()

d_df_ref = pd.DataFrame(d_ref, columns = columns_d.keys()).astype(columns_d).set_index('parcel_id')
d_df_ref['geometry'] = d_df_ref['geometry_wkt'].apply(loads)
d_df_ref = gpd.GeoDataFrame(d_df_ref, geometry='geometry', crs="EPSG:28992")  # Assuming WGS84 (latitude/longitude)
d_df_ref.drop(columns=['geometry_wkt'], inplace=True)


sql_select_query_ref = """SELECT perceel_id, dekking_ve,dekking_mo, peil_gemid, peil_gem_1 FROM a_input.a_2016"""

cursor.execute(sql_select_query_ref)
a_ref = cursor.fetchall()

a_df_ref = pd.DataFrame(a_ref, columns = columns_a.keys()).astype(columns_a).set_index('perceel_id')

#%%

d_df['dekking_ve'] = d_df.index.map(a_df['dekking_ve'])
d_df['dekking_mo'] = d_df.index.map(a_df['dekking_mo'])
  
d_df['dekking_org_totaal'] = (d_df['dekking_ve'] + d_df['dekking_mo']) / 100
d_df['median_per_ha'] = d_df['median'] / ((d_df.geometry.area / 10000) * d_df['dekking_org_totaal'])

d_df_ref['dekking_ve'] = d_df_ref.index.map(a_df_ref['dekking_ve'])
d_df_ref['dekking_mo'] = d_df_ref.index.map(a_df_ref['dekking_mo'])
  
d_df_ref['dekking_org_totaal'] = (d_df_ref['dekking_ve'] + d_df_ref['dekking_mo']) / 100
d_df_ref['median_per_ha'] = d_df_ref['median'] / ((d_df_ref.geometry.area / 10000) * d_df_ref['dekking_org_totaal'])

d_df['org_areaal'] = d_df['dekking_org_totaal'] * d_df.geometry.area

d_df_ref['org_areaal'] = d_df_ref['dekking_org_totaal'] * d_df_ref.geometry.area #m2

#%%
provincies = gpd.read_file(r'P:/11207812-somers-ontwikkeling/base_data/GIS_layers/Bestuurlijke_grenzen/B1_Provinciale_indeling_van_NederlandPolygon.shp')

provincie_ls = []

for parcel_id, row in d_df_ref.iterrows():
    #print (row)
    geom_row = row.geometry.centroid
    province = provincies[provincies.contains(geom_row)]['PROVC_NM']
    if province.empty:
        province = np.nan
    else:
        province = province.iloc[0]
    provincie_ls.append(province)

d_df_ref['provincie'] = provincie_ls

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


#%% Merge peildata
a_df = a_df.add_suffix('_2023')
a_df_ref = a_df_ref.add_suffix('_2016')

merged_df = a_df.merge(a_df_ref[['peil_gemid_2016', 'peil_gem_1_2016']], left_index=True, right_index=True, how='left')
merged_df = merged_df[merged_df.index.isin(d_df.index)]

#%% Select infiltration measures

infiltratie_mask = d_df['measure'].isin(['SSI', 'PSSI'])
infiltratie_percelen = d_df[infiltratie_mask]

#%% Select peilverhoging, minstens 5cm

peilverhoging_percelen = merged_df[
    (merged_df['peil_gemid_2023'] >= 0.05 + merged_df['peil_gemid_2016']) |
    (merged_df['peil_gem_1_2023'] >= 0.05 + merged_df['peil_gem_1_2016'])]

#%% Select peilverlaging, minstens 5cm 
 
peilverlaging_percelen = merged_df[
    (merged_df['peil_gemid_2023'] <= -0.05 + merged_df['peil_gemid_2016'] ) |
    (merged_df['peil_gem_1_2023'] <= -0.05 + merged_df['peil_gem_1_2016'] )]

#%% Uitstoot maatregel percelen

infiltratie_percelen_ref = d_df_ref.loc[d_df_ref.index.isin(infiltratie_percelen.index)]
peilverlaging_percelen_ref =  d_df_ref.loc[d_df_ref.index.isin(peilverlaging_percelen.index)]
peilverhoging_percelen_ref =  d_df_ref.loc[d_df_ref.index.isin(peilverhoging_percelen.index)]

peilverlaging_percelen_monitoringsjaar =  d_df.loc[d_df.index.isin(peilverlaging_percelen.index)]
peilverhoging_percelen_monitoringsjaar =  d_df.loc[d_df.index.isin(peilverhoging_percelen.index)]

print("Infiltratie percelen 2016", f"{infiltratie_percelen_ref['median'].sum():.0f}")
print(f"Infiltratie percelen {monitoringsjaar}",f"{infiltratie_percelen['median'].sum():.0f}")
print("Peilverlaging percelen 2016",f"{peilverlaging_percelen_ref['median'].sum():.0f}")
print(f"Peilverlaging percelen {monitoringsjaar}",f"{peilverlaging_percelen_monitoringsjaar['median'].sum():.0f}")
print("Peilverhoging percelen 2016",f"{peilverhoging_percelen_ref['median'].sum():.0f}")
print(f"Peilverhoging percelen {monitoringsjaar}",f"{peilverhoging_percelen_monitoringsjaar['median'].sum():.0f}")
#%% Arealen maatregelen

area_inf = infiltratie_percelen['org_areaal'].sum () / 10000 #ha
area_peilverhoging = peilverhoging_percelen_ref['org_areaal'].sum () / 10000 #ha
area_peilverlaging = peilverlaging_percelen_ref['org_areaal'].sum () / 10000 #ha

print(f"Areaal infiltratie percelen {monitoringsjaar}", f"{area_inf:.2f}")
print(f"Areaal peilverhoging percelen {monitoringsjaar}",f"{area_peilverhoging:.2f}")
print(f"Areaal peilverlaging percelen {monitoringsjaar}",f"{area_peilverlaging:.2f}")

#%% Infiltratie per provincie monitoringsjaar

inf_zuid_holland = infiltratie_percelen[infiltratie_percelen['provincie'] == 'Zuid-Holland']
inf_utrecht = infiltratie_percelen[infiltratie_percelen['provincie'] == 'Utrecht']
inf_noord_holland = infiltratie_percelen[infiltratie_percelen['provincie'] == 'Noord-Holland']
inf_overijssel = infiltratie_percelen[infiltratie_percelen['provincie'] == 'Overijssel']
inf_drenthe = infiltratie_percelen[infiltratie_percelen['provincie'] == 'Drenthe']
inf_groningen = infiltratie_percelen[infiltratie_percelen['provincie'] == 'Groningen']
inf_friesland = infiltratie_percelen[infiltratie_percelen['provincie'] == 'Friesland']
inf_gelderland = infiltratie_percelen[infiltratie_percelen['provincie'] == 'Gelderland']
inf_zeeland = infiltratie_percelen[infiltratie_percelen['provincie'] == 'Zeeland']
inf_flevoland = infiltratie_percelen[infiltratie_percelen['provincie'] == 'Flevoland']
inf_brabant = infiltratie_percelen[infiltratie_percelen['provincie'] == 'Noord-Brabant']

print('Total emissisions in groningen:', f"{inf_groningen['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in friesland:', f"{inf_friesland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in overijssel:', f"{inf_overijssel['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in drenthe:', f"{inf_drenthe['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in gelderland:', f"{inf_gelderland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in utrecht:', f"{inf_utrecht['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in noorinf_holland:', f"{inf_noord_holland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in zeeland:', f"{inf_zeeland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in zuiinf_holland:', f"{inf_zuid_holland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in flevoland:', f"{inf_flevoland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in brabant:', f"{inf_brabant['median'].sum() / (1000):.0f}", 'ton CO2')

print('Total areaal inf groningen:', f"{inf_groningen.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal inf friesland:', f"{inf_friesland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal inf overijssel:', f"{inf_overijssel.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal inf drenthe:', f"{inf_drenthe.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal inf gelderland:', f"{inf_gelderland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal inf utrecht:', f"{inf_utrecht.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal inf noorinf_holland:', f"{inf_noord_holland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal inf zeeland:', f"{inf_zeeland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal inf zuiinf_holland:', f"{inf_zuid_holland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal inf flevoland:', f"{inf_flevoland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal inf brabant:', f"{inf_brabant.org_areaal.sum() / 10000:.0f}", 'ha')


#%% Peilverhoging per provincie monitoringsjaar

verhoging_zuid_holland = peilverhoging_percelen_monitoringsjaar[peilverhoging_percelen_monitoringsjaar['provincie'] == 'Zuid-Holland']
verhoging_utrecht = peilverhoging_percelen_monitoringsjaar[peilverhoging_percelen_monitoringsjaar['provincie'] == 'Utrecht']
verhoging_noord_holland = peilverhoging_percelen_monitoringsjaar[peilverhoging_percelen_monitoringsjaar['provincie'] == 'Noord-Holland']
verhoging_overijssel = peilverhoging_percelen_monitoringsjaar[peilverhoging_percelen_monitoringsjaar['provincie'] == 'Overijssel']
verhoging_drenthe = peilverhoging_percelen_monitoringsjaar[peilverhoging_percelen_monitoringsjaar['provincie'] == 'Drenthe']
verhoging_groningen = peilverhoging_percelen_monitoringsjaar[peilverhoging_percelen_monitoringsjaar['provincie'] == 'Groningen']
verhoging_friesland = peilverhoging_percelen_monitoringsjaar[peilverhoging_percelen_monitoringsjaar['provincie'] == 'Friesland']
verhoging_gelderland = peilverhoging_percelen_monitoringsjaar[peilverhoging_percelen_monitoringsjaar['provincie'] == 'Gelderland']
verhoging_zeeland = peilverhoging_percelen_monitoringsjaar[peilverhoging_percelen_monitoringsjaar['provincie'] == 'Zeeland']
verhoging_flevoland = peilverhoging_percelen_monitoringsjaar[peilverhoging_percelen_monitoringsjaar['provincie'] == 'Flevoland']
verhoging_brabant = peilverhoging_percelen_monitoringsjaar[peilverhoging_percelen_monitoringsjaar['provincie'] == 'Noord-Brabant']

print('Total emissisions in groningen:', f"{verhoging_groningen['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in friesland:', f"{verhoging_friesland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in overijssel:', f"{verhoging_overijssel['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in drenthe:', f"{verhoging_drenthe['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in gelderland:', f"{verhoging_gelderland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in utrecht:', f"{verhoging_utrecht['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in noord_holland:', f"{verhoging_noord_holland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in zeeland:', f"{verhoging_zeeland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in zuid_holland:', f"{verhoging_zuid_holland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in flevoland:', f"{verhoging_flevoland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in brabant:', f"{verhoging_brabant['median'].sum() / (1000):.0f}", 'ton CO2')

print('Total areaal peilverhoging groningen:', f"{verhoging_groningen.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverhoging friesland:', f"{verhoging_friesland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverhoging overijssel:', f"{verhoging_overijssel.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverhoging drenthe:', f"{verhoging_drenthe.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverhoging gelderland:', f"{verhoging_gelderland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverhoging utrecht:', f"{verhoging_utrecht.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverhoging noord_holland:', f"{verhoging_noord_holland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverhoging zeeland:', f"{verhoging_zeeland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverhoging zuid_holland:', f"{verhoging_zuid_holland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverhoging flevoland:', f"{verhoging_flevoland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverhoging brabant:', f"{verhoging_brabant.org_areaal.sum() / 10000:.0f}", 'ha')

#%% Peilverlaging per provincie monitoringsjaar

verlaging_zuid_holland = peilverlaging_percelen_monitoringsjaar[peilverlaging_percelen_monitoringsjaar['provincie'] == 'Zuid-Holland']
verlaging_utrecht = peilverlaging_percelen_monitoringsjaar[peilverlaging_percelen_monitoringsjaar['provincie'] == 'Utrecht']
verlaging_noord_holland = peilverlaging_percelen_monitoringsjaar[peilverlaging_percelen_monitoringsjaar['provincie'] == 'Noord-Holland']
verlaging_overijssel = peilverlaging_percelen_monitoringsjaar[peilverlaging_percelen_monitoringsjaar['provincie'] == 'Overijssel']
verlaging_drenthe = peilverlaging_percelen_monitoringsjaar[peilverlaging_percelen_monitoringsjaar['provincie'] == 'Drenthe']
verlaging_groningen = peilverlaging_percelen_monitoringsjaar[peilverlaging_percelen_monitoringsjaar['provincie'] == 'Groningen']
verlaging_friesland = peilverlaging_percelen_monitoringsjaar[peilverlaging_percelen_monitoringsjaar['provincie'] == 'Friesland']
verlaging_gelderland = peilverlaging_percelen_monitoringsjaar[peilverlaging_percelen_monitoringsjaar['provincie'] == 'Gelderland']
verlaging_zeeland = peilverlaging_percelen_monitoringsjaar[peilverlaging_percelen_monitoringsjaar['provincie'] == 'Zeeland']
verlaging_flevoland = peilverlaging_percelen_monitoringsjaar[peilverlaging_percelen_monitoringsjaar['provincie'] == 'Flevoland']
verlaging_brabant = peilverlaging_percelen_monitoringsjaar[peilverlaging_percelen_monitoringsjaar['provincie'] == 'Noord-Brabant']

print('Total emissisions in groningen:', f"{verlaging_groningen['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in friesland:', f"{verlaging_friesland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in overijssel:', f"{verlaging_overijssel['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in drenthe:', f"{verlaging_drenthe['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in gelderland:', f"{verlaging_gelderland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in utrecht:', f"{verlaging_utrecht['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in noord_holland:', f"{verlaging_noord_holland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in zeeland:', f"{verlaging_zeeland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in zuid_holland:', f"{verlaging_zuid_holland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in flevoland:', f"{verlaging_flevoland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in brabant:', f"{verlaging_brabant['median'].sum() / (1000):.0f}", 'ton CO2')

print('Total areaal peilverlaging groningen:', f"{verlaging_groningen.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverlaging friesland:', f"{verlaging_friesland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverlaging overijssel:', f"{verlaging_overijssel.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverlaging drenthe:', f"{verlaging_drenthe.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverlaging gelderland:', f"{verlaging_gelderland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverlaging utrecht:', f"{verlaging_utrecht.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverlaging noord_holland:', f"{verlaging_noord_holland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverlaging zeeland:', f"{verlaging_zeeland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverlaging zuid_holland:', f"{verlaging_zuid_holland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverlaging flevoland:', f"{verlaging_flevoland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverlaging brabant:', f"{verlaging_brabant.org_areaal.sum() / 10000:.0f}", 'ha')

#%% Infiltratie per provincie refjaar

inf_zuid_holland = infiltratie_percelen_ref[infiltratie_percelen_ref['provincie'] == 'Zuid-Holland']
inf_utrecht = infiltratie_percelen_ref[infiltratie_percelen_ref['provincie'] == 'Utrecht']
inf_noord_holland = infiltratie_percelen_ref[infiltratie_percelen_ref['provincie'] == 'Noord-Holland']
inf_overijssel = infiltratie_percelen_ref[infiltratie_percelen_ref['provincie'] == 'Overijssel']
inf_drenthe = infiltratie_percelen_ref[infiltratie_percelen_ref['provincie'] == 'Drenthe']
inf_groningen = infiltratie_percelen_ref[infiltratie_percelen_ref['provincie'] == 'Groningen']
inf_friesland = infiltratie_percelen_ref[infiltratie_percelen_ref['provincie'] == 'Friesland']
inf_gelderland = infiltratie_percelen_ref[infiltratie_percelen_ref['provincie'] == 'Gelderland']
inf_zeeland = infiltratie_percelen_ref[infiltratie_percelen_ref['provincie'] == 'Zeeland']
inf_flevoland = infiltratie_percelen_ref[infiltratie_percelen_ref['provincie'] == 'Flevoland']
inf_brabant = infiltratie_percelen_ref[infiltratie_percelen_ref['provincie'] == 'Noord-Brabant']

print('Total emissisions in groningen:', f"{inf_groningen['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in friesland:', f"{inf_friesland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in overijssel:', f"{inf_overijssel['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in drenthe:', f"{inf_drenthe['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in gelderland:', f"{inf_gelderland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in utrecht:', f"{inf_utrecht['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in noorinf_holland:', f"{inf_noord_holland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in zeeland:', f"{inf_zeeland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in zuiinf_holland:', f"{inf_zuid_holland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in flevoland:', f"{inf_flevoland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in brabant:', f"{inf_brabant['median'].sum() / (1000):.0f}", 'ton CO2')

print('Total areaal inf groningen:', f"{inf_groningen.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal inf friesland:', f"{inf_friesland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal inf overijssel:', f"{inf_overijssel.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal inf drenthe:', f"{inf_drenthe.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal inf gelderland:', f"{inf_gelderland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal inf utrecht:', f"{inf_utrecht.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal inf noorinf_holland:', f"{inf_noord_holland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal inf zeeland:', f"{inf_zeeland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal inf zuiinf_holland:', f"{inf_zuid_holland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal inf flevoland:', f"{inf_flevoland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal inf brabant:', f"{inf_brabant.org_areaal.sum() / 10000:.0f}", 'ha')


#%% Peilverhoging per provincie refjaar

verhoging_zuid_holland = peilverhoging_percelen_ref[peilverhoging_percelen_ref['provincie'] == 'Zuid-Holland']
verhoging_utrecht = peilverhoging_percelen_ref[peilverhoging_percelen_ref['provincie'] == 'Utrecht']
verhoging_noord_holland = peilverhoging_percelen_ref[peilverhoging_percelen_ref['provincie'] == 'Noord-Holland']
verhoging_overijssel = peilverhoging_percelen_ref[peilverhoging_percelen_ref['provincie'] == 'Overijssel']
verhoging_drenthe = peilverhoging_percelen_ref[peilverhoging_percelen_ref['provincie'] == 'Drenthe']
verhoging_groningen = peilverhoging_percelen_ref[peilverhoging_percelen_ref['provincie'] == 'Groningen']
verhoging_friesland = peilverhoging_percelen_ref[peilverhoging_percelen_ref['provincie'] == 'Friesland']
verhoging_gelderland = peilverhoging_percelen_ref[peilverhoging_percelen_ref['provincie'] == 'Gelderland']
verhoging_zeeland = peilverhoging_percelen_ref[peilverhoging_percelen_ref['provincie'] == 'Zeeland']
verhoging_flevoland = peilverhoging_percelen_ref[peilverhoging_percelen_ref['provincie'] == 'Flevoland']
verhoging_brabant = peilverhoging_percelen_ref[peilverhoging_percelen_ref['provincie'] == 'Noord-Brabant']

print('Total emissisions in groningen:', f"{verhoging_groningen['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in friesland:', f"{verhoging_friesland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in overijssel:', f"{verhoging_overijssel['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in drenthe:', f"{verhoging_drenthe['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in gelderland:', f"{verhoging_gelderland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in utrecht:', f"{verhoging_utrecht['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in noord_holland:', f"{verhoging_noord_holland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in zeeland:', f"{verhoging_zeeland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in zuid_holland:', f"{verhoging_zuid_holland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in flevoland:', f"{verhoging_flevoland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in brabant:', f"{verhoging_brabant['median'].sum() / (1000):.0f}", 'ton CO2')

print('Total areaal peilverhoging groningen:', f"{verhoging_groningen.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverhoging friesland:', f"{verhoging_friesland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverhoging overijssel:', f"{verhoging_overijssel.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverhoging drenthe:', f"{verhoging_drenthe.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverhoging gelderland:', f"{verhoging_gelderland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverhoging utrecht:', f"{verhoging_utrecht.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverhoging noord_holland:', f"{verhoging_noord_holland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverhoging zeeland:', f"{verhoging_zeeland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverhoging zuid_holland:', f"{verhoging_zuid_holland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverhoging flevoland:', f"{verhoging_flevoland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverhoging brabant:', f"{verhoging_brabant.org_areaal.sum() / 10000:.0f}", 'ha')

#%% Peilverlaging per provincie refjaar

verlaging_zuid_holland = peilverlaging_percelen_ref[peilverlaging_percelen_ref['provincie'] == 'Zuid-Holland']
verlaging_utrecht = peilverlaging_percelen_ref[peilverlaging_percelen_ref['provincie'] == 'Utrecht']
verlaging_noord_holland = peilverlaging_percelen_ref[peilverlaging_percelen_ref['provincie'] == 'Noord-Holland']
verlaging_overijssel = peilverlaging_percelen_ref[peilverlaging_percelen_ref['provincie'] == 'Overijssel']
verlaging_drenthe = peilverlaging_percelen_ref[peilverlaging_percelen_ref['provincie'] == 'Drenthe']
verlaging_groningen = peilverlaging_percelen_ref[peilverlaging_percelen_ref['provincie'] == 'Groningen']
verlaging_friesland = peilverlaging_percelen_ref[peilverlaging_percelen_ref['provincie'] == 'Friesland']
verlaging_gelderland = peilverlaging_percelen_ref[peilverlaging_percelen_ref['provincie'] == 'Gelderland']
verlaging_zeeland = peilverlaging_percelen_ref[peilverlaging_percelen_ref['provincie'] == 'Zeeland']
verlaging_flevoland = peilverlaging_percelen_ref[peilverlaging_percelen_ref['provincie'] == 'Flevoland']
verlaging_brabant = peilverlaging_percelen_ref[peilverlaging_percelen_ref['provincie'] == 'Noord-Brabant']

print('Total emissisions in groningen:', f"{verlaging_groningen['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in friesland:', f"{verlaging_friesland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in overijssel:', f"{verlaging_overijssel['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in drenthe:', f"{verlaging_drenthe['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in gelderland:', f"{verlaging_gelderland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in utrecht:', f"{verlaging_utrecht['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in noord_holland:', f"{verlaging_noord_holland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in zeeland:', f"{verlaging_zeeland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in zuid_holland:', f"{verlaging_zuid_holland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in flevoland:', f"{verlaging_flevoland['median'].sum() / (1000):.0f}", 'ton CO2')
print('Total emissisions in brabant:', f"{verlaging_brabant['median'].sum() / (1000):.0f}", 'ton CO2')

print('Total areaal peilverlaging groningen:', f"{verlaging_groningen.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverlaging friesland:', f"{verlaging_friesland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverlaging overijssel:', f"{verlaging_overijssel.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverlaging drenthe:', f"{verlaging_drenthe.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverlaging gelderland:', f"{verlaging_gelderland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverlaging utrecht:', f"{verlaging_utrecht.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverlaging noord_holland:', f"{verlaging_noord_holland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverlaging zeeland:', f"{verlaging_zeeland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverlaging zuid_holland:', f"{verlaging_zuid_holland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverlaging flevoland:', f"{verlaging_flevoland.org_areaal.sum() / 10000:.0f}", 'ha')
print('Total areaal peilverlaging brabant:', f"{verlaging_brabant.org_areaal.sum() / 10000:.0f}", 'ha')

#%% Check adminsitratieve veranderingen: alle veranderingen en missende percelen 

peilverhoging_percelen = merged_df[
    (merged_df['peil_gemid_2023'] > merged_df['peil_gemid_2016'])|
    (merged_df['peil_gem_1_2023'] > merged_df['peil_gem_1_2016']) 
    & ~merged_df.index.isin(infiltratie_percelen.index)]

peilverlaging_percelen = merged_df[
    (merged_df['peil_gemid_2023'] < merged_df['peil_gemid_2016'] ) |
    (merged_df['peil_gem_1_2023'] < merged_df['peil_gem_1_2016'] ) 
    & ~merged_df.index.isin(infiltratie_percelen.index)
    & ~merged_df.index.isin(peilverhoging_percelen.index)]

peilverlaging_percelen_ref = d_df_ref.loc[d_df_ref.index.isin(peilverlaging_percelen.index)]
peilverhoging_percelen_ref = d_df_ref.loc[d_df_ref.index.isin(peilverhoging_percelen.index)]

peilverlaging_percelen_monitoringsjaar = d_df.loc[d_df.index.isin(peilverlaging_percelen.index)]
peilverhoging_percelen_monitoringsjaar = d_df.loc[d_df.index.isin(peilverhoging_percelen.index)]

missende_percelen = d_df_ref.loc[~d_df_ref.index.isin(d_df.index)]
extra_percelen = d_df.loc[~d_df.index.isin(d_df_ref.index)]


print("Missende percelen 2023",f"{missende_percelen['median'].sum():.0f}")
print("Extra percelen 2023",f"{extra_percelen['median'].sum():.0f}")
print("Peilverlaging percelen 2016",f"{peilverlaging_percelen_ref['median'].sum():.0f}")
print(f"Peilverlaging percelen {monitoringsjaar}",f"{peilverlaging_percelen_monitoringsjaar['median'].sum():.0f}")
print("Peilverhoging percelen 2016",f"{peilverhoging_percelen_ref['median'].sum():.0f}")
print(f"Peilverhoging percelen {monitoringsjaar}",f"{peilverhoging_percelen_monitoringsjaar['median'].sum():.0f}")