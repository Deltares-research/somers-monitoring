import configparser
import psycopg2
import pandas as pd
from scipy.stats import gaussian_kde
import numpy as np
from shapely import wkb, wkt
from shapely.wkt import loads
import geopandas as gpd
from sklearn.linear_model import LinearRegression

import somers
from somers import to_linux
import os

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors

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
#%%

sql_select_query = f"""SELECT parcel_id, minimum, median, maximum, perc_10, 
    perc_90, measure, geometry, ST_AsText(geometry) FROM d_post.d_{monitoringsjaar}"""

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
            'geometry_wkt':object,}

d_df = pd.DataFrame(d_monitoringsjaar, columns = columns_d.keys()).astype(columns_d).set_index('parcel_id')
d_df['geometry'] = d_df['geometry_wkt'].apply(loads)
d_df = gpd.GeoDataFrame(d_df, geometry='geometry', crs="EPSG:28992") 
d_df.drop(columns=['geometry_wkt'], inplace=True)



sql_select_query = f"""SELECT perceel_id, peil_gemid, ahn3_media, ahn5_media, perceelsbr, archetype_, dekking_ve,dekking_mo FROM a_input.a_{monitoringsjaar}"""

cursor.execute(sql_select_query)
a_monitoringsjaar = cursor.fetchall()

columns_a = {'perceel_id': str,
             'peil_gemid': float,
             'ahn3_media': float,
             'ahn5_media': float,
             'perceelsbr': float,
             'archetype_': str,
             'dekking_ve': float,
             'dekking_mo': float,}

a_df = pd.DataFrame(a_monitoringsjaar, columns = columns_a.keys()).astype(columns_a).set_index('perceel_id')


sql_select_query = f"""SELECT parcel_id, mean_summer_gwl, mean_year_gwl FROM c_output.c_{monitoringsjaar}"""

cursor.execute(sql_select_query)
c_monitoringsjaar = cursor.fetchall()

columns_c = {'parcel_id': str,
             'mean_summer_gwl': float,
             'mean_year_gwl': float,}

c_df = pd.DataFrame(c_monitoringsjaar, columns = columns_c.keys()).astype(columns_c).set_index('parcel_id')

#%% 

#add freeboard and ditch distance
d_df = d_df.merge(a_df[['peil_gemid']], left_index=True, right_index=True, how='left')
d_df = d_df.merge(a_df[['ahn3_media']], left_index=True, right_index=True, how='left') #AHN 3!
d_df = d_df.merge(a_df[['perceelsbr']], left_index=True, right_index=True, how='left')
d_df = d_df.merge(a_df[['archetype_']], left_index=True, right_index=True, how='left')
d_df = d_df.merge(a_df[['dekking_ve']], left_index=True, right_index=True, how='left')
d_df = d_df.merge(a_df[['dekking_mo']], left_index=True, right_index=True, how='left')
d_df = d_df.merge(c_df[['mean_summer_gwl']], left_index=True, right_index=True, how='left') #m-mv
d_df = d_df.merge(c_df[['mean_year_gwl']], left_index=True, right_index=True, how='left') #m-mv
d_df['dekkingsgraad_totaal'] = (d_df['dekking_ve'] + d_df['dekking_mo']) / 100
d_df.rename(columns = {'peil_gemid':'summer_stage', 'ahn3_media':'surface_level', 'perceelsbr':'parcel_width'}, inplace=True)
d_df['freeboard'] = d_df['surface_level'] - d_df['summer_stage']
d_df['median_per_ha'] = d_df['median'] / ((d_df.geometry.area / 10000) * d_df['dekkingsgraad_totaal'])

#%% Monitoringsrapport plot
# plot t.o.v. gws en perceelsbreedte

plt.figure(figsize=(8, 6))
scatter_1 = plt.scatter(d_df['mean_summer_gwl'], d_df['median_per_ha'], c=d_df['parcel_width'], cmap='viridis')
cbar_1 = plt.colorbar(scatter_1)
cbar_1.set_label('Perceelsbreedte [m]')
plt.xlabel('Gem. zomer grondwaterstand [m-mv]')
plt.ylabel('Mediane uitstoot [kg CO2/ha/jaar]')
plt.title('Berekende uitstoot uit kustvlakte percelen')
plt.xlim(0,2.5)

#%% Plot t.o.v. gws en archetype

soil_map = {
    'Vc': 'hVb', 'aVc': 'hVb', 'hVc': 'hVb', 'hVb': 'hVb',
    'pVb': 'pVb', 'pVc': 'pVb', 'kVc': 'pVb',
    'Vk': 'hVk', 'hVk': 'hVk',
    'hVz': 'hVz', 'aVz': 'hVz', 'Vz': 'hVz',
    'Vp': 'Vp', 'iVp': 'Vp',
    'Wo': 'Wo', 'hVs': 'hVs', 'AVk':'AVk',
    'iWp': 'iWp', 'kVk': 'kVk', 'kVs': 'kVs', 'kVz': 'kVz',
    'kWp': 'kWp', 'vWz': 'vWz', 'zVp': 'zVp', 'zVz': 'zVz', 'zWp': 'zWp'}

custom_order = ['kVk', 'pVb', 'kVs', 'kVz', 'hVk', 'hVb', 'hVs', 'hVz', 'Vp', 'zVz', 'zVp', 'AVk', 'kWp', 'Wo', 'vWz', 'iWp', 'zWp']
d_df['soil_group'] = pd.Categorical(d_df['archetype_'], categories=custom_order, ordered=True)

d_df['soil_group'] = d_df['soil_group'].astype('category')
soil_labels = d_df['soil_group'].cat.categories
soil_codes = d_df['soil_group'].cat.codes

cmap = plt.cm.get_cmap('tab20b', 17)
norm = mcolors.BoundaryNorm(np.arange(-0.5, 17.5, 1), cmap.N)

plt.figure(figsize=(10, 7))
sc = plt.scatter(
    d_df['mean_summer_gwl'],
    d_df['median_per_ha'],
    c=soil_codes,
    cmap=cmap,
    norm=norm,
    s=50
)

cbar = plt.colorbar(sc, ticks=np.arange(17))
cbar.ax.set_yticklabels(soil_labels)
cbar.set_label('Subklasse')

plt.xlabel('Gem. zomer grondwaterstand [m-mv]')
plt.ylabel('Mediane uitstoot [kg CO2/ha/jaar]')
plt.title('Berekende uitstoot uit kustvlakte percelen')
plt.xlim(0, 2.5)
plt.tight_layout()
plt.show()

#%% Plot t.o.v. provincies

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


# Map province names to integer codes for coloring
prov_labels = sorted(d_df['provincie'].dropna().unique())
prov_map = {prov: i for i, prov in enumerate(prov_labels)}
d_df['prov_code'] = d_df['provincie'].map(prov_map)

# Plotting
cmap = plt.cm.get_cmap('tab20b', len(prov_labels))
norm = mcolors.BoundaryNorm(np.arange(-0.5, len(prov_labels) + 0.5, 1), cmap.N)

plt.figure(figsize=(12, 7))
sc = plt.scatter(
    d_df['mean_summer_gwl'],
    d_df['median_per_ha'],
    c=d_df['prov_code'],
    cmap=cmap,
    norm=norm,
    s=50
)

cbar = plt.colorbar(sc, ticks=np.arange(len(prov_labels)))
cbar.ax.set_yticklabels(prov_labels)
cbar.set_label('Provincie')

plt.xlabel('Gem. zomer grondwaterstand [m-mv]')
plt.ylabel('Mediane uitstoot [kg CO2/ha/jaar]')
plt.title('Berekende uitstoot uit kustvlakte percelen')
plt.xlim(0, 2.5)
plt.tight_layout()
plt.show()


#%% Plot t.o.v. gws, density plot

x = d_df['mean_summer_gwl']
y = d_df['median_per_ha']

mask = np.isfinite(x) & np.isfinite(y)
x_clean = x[mask]
y_clean = y[mask]

xy = np.vstack([x_clean, y_clean])
z = gaussian_kde(xy)(xy)

fig, ax = plt.subplots(figsize=(8, 6))
sc = ax.scatter(x_clean, y_clean, c=z, s=50, cmap='viridis')
plt.colorbar(sc, label='Puntdichtheid')
ax.set_xlabel('Gem. zomer grondwaterstand [m-mv]')
ax.set_ylabel('Mediane uitstoot [kg CO2/ha/jaar]')
ax.set_title('Berekende uitstoot uit kustvlakte percelen')
ax.set_xlim(0, 2.5)
plt.show()


#%% Plot high emission and high groundwaterlevel

df = d_df[['mean_summer_gwl', 'median_per_ha']].dropna()

X = df['mean_summer_gwl'].values.reshape(-1, 1)
y = df['median_per_ha'].values
model = LinearRegression()
model.fit(X, y)
df['trend_pred'] = model.predict(X)
df['residual'] = df['median_per_ha'] - df['trend_pred']

threshold = df['residual'].quantile(0.98)
top_10_df = df[df['residual'] >= threshold]

filtered_df = top_10_df[top_10_df['mean_summer_gwl'] <= 0.81]

filtered_full = d_df.loc[filtered_df.index]

plt.figure(figsize=(10, 6))
plt.scatter(df['mean_summer_gwl'], df['median_per_ha'], color='lightgray', s=10, label='All data')
plt.scatter(filtered_df['mean_summer_gwl'], filtered_df['median_per_ha'], color='red', s=15, label='Top 10% & GWL ≤ 1.0m')

x_vals = np.linspace(0, 2.5, 100).reshape(-1, 1)
y_vals = model.predict(x_vals)
plt.plot(x_vals, y_vals, 'k--', label='Trend line')

plt.xlabel('Gem. zomer grondwaterstand [m-mv]')
plt.ylabel('Mediane uitstoot [kg CO2/ha/jaar]')
plt.title('Top 10% Above Trend with GWL ≤ 1.0m')
plt.legend()
plt.tight_layout()
plt.show()

#filtered_full.to_file(r'P:/11207812-somers-uitvoering/monitoring_2023/src/analyse_src/hoge_uitstoot_hoge_gws.gpkg')
