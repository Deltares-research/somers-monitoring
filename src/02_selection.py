import pandas as pd
import warnings
import numpy as np
import configparser
import psycopg2
import sys
import shapely
import os

os.environ["USE_PYGEOS"] = "0"
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
warnings.filterwarnings("ignore", category=DeprecationWarning)

from somers import to_linux

sys.path.append(to_linux(r"P:/11207812-somers-uitvoering/monitoring_2023/src/"))

# %%
# Establish database connection & make table
user_credentials = "ln"
config_file_path = to_linux(
    f"P:/11207812-somers-uitvoering/monitoring_2023/configfile_somers_db_{user_credentials}.txt"
)

# %%

# Read the configuration file
cf = configparser.ConfigParser()
cf.read(config_file_path)

# Get the database connection parameters from the configuration file
host = cf.get("PostGIS", "host")
database = cf.get("PostGIS", "db")
user = cf.get("PostGIS", "user")
password = cf.get("PostGIS", "pass")
port = cf.get("PostGIS", "port")

# Establish a connection to the PostgreSQL database
connection = psycopg2.connect(
    host=host, database=database, user=user, password=password, port=port
)
cursor = connection.cursor()
# %%

monitoringsjaar = 2022
# %%

sql_create_table = f"""
    CREATE TABLE IF NOT EXISTS b_selection.b_{monitoringsjaar}(
        parcel_id TEXT UNIQUE NOT NULL,
        x FLOAT NOT NULL,
        y FLOAT NOT NULL,
        archetype TEXT NOT NULL,
        surface_level_ahn3 FLOAT NOT NULL,
        surface_level_ahn5 FLOAT NOT NULL,
        parcel_width FLOAT NOT NULL,
        summer_stage FLOAT NOT NULL,
        winter_stage FLOAT NOT NULL,
        SSI_distance FLOAT,
        SSI_depth FLOAT,
        PSSI_distance FLOAT,
        PSSI_depth FLOAT,
        measure TEXT,
        geometry GEOMETRY);
    
        ALTER TABLE b_selection.b_{monitoringsjaar} OWNER TO somers
"""

cursor.execute(sql_create_table)
connection.commit()

# Download a
sql_query = f"""SELECT *, ST_AsText(geom) FROM a_input.a_{monitoringsjaar}  ORDER BY perceel_id"""
cursor.execute(sql_query)
records = cursor.fetchall()

columns = {
    "gid": int,
    "perceel_id": str,
    "type": str,  # BGT perceelstype
    "dekking_ve": int,  # Dekking veen
    "dekking_mo": int,  # dekking moerig
    "dekking_be": int,  # dekking begraven veen
    "archetype_": str,  # hoofd archetype
    "archetyp_1": str,  # tweede archetype
    "meststoffe": str,  # wel of niet meststoffenperceel
    "opmerking_": str,
    "perceelsbr": float,  # perceelsbreedte
    "ahn3_media": float,  # ahn3
    "ahn4_media": float,  # ahn4
    "ahn5_media": float,  # ahn5
    "peil_gemid": float,  # zomerpeil
    "peil_gem_1": float,  # winterpeil
    "peil_aanta": float,  # peilaantal
    "peil_gebie": float,  # Peil_Gebied1_zp
    "peil_geb_1": float,  # Peil_Gebied1_wp
    "peil_geb_2": float,  # Peil_Gebied2_zp
    "peil_geb_3": float,  # Peil_Gebied2_wp
    "maatregel": str,  # wel of geen maatregel
    "owd_draina": float,  # OWD drainafstand
    "owd_draind": float,  # OWD draindiepte
    "dd_scenari": str,  # scenario druk drainage
    "dd_drainaf": float,  # DD drainafstand
    "dd_draindi": float,  # DD draindiepte
    "drain_perc": float,  # percentage perceel gedekt door maatregel
    "geometry": object,
    "geometry_wkt": object,
}

dataframe = (
    pd.DataFrame(records, columns=columns.keys())
    .astype(columns)
    .set_index("perceel_id")
)
dataframe["geometry"] = dataframe["geometry_wkt"].apply(shapely.wkt.loads)

# %%


def measure_func(SSI_dist, PSSI_dist):
    if not np.isnan(SSI_dist):
        measure = "SSI"
    elif not np.isnan(PSSI_dist):
        measure = "PSSI"
    else:
        measure = "ref"
    return measure


for n, row in dataframe.iterrows():
    print(n)
    # Unpack row
    parcel_id = str(n)
    x = float(row.geometry.centroid.x)
    y = float(row.geometry.centroid.y)
    archetype = str(row.archetype_)
    bgt_type = str(row.type)
    surface_level_ahn3 = float(row.ahn3_media) if row.ahn3_media is not None else None
    surface_level_ahn5 = (
        float(row.ahn5_media) if row.ahn5_media != 8888 else row.ahn3_media
    )
    parcel_width = float(row.perceelsbr) if row.perceelsbr is not None else None
    summer_stage = float(row.peil_gemid) if row.peil_gemid is not None else None
    winter_stage = float(row.peil_gem_1) if row.peil_gem_1 is not None else None
    SSI_distance = float(row.owd_draina) if row.owd_draina is not None else None
    SSI_depth = float(row.owd_draind) if row.owd_draind is not None else None
    PSSI_distance = float(row.dd_drainaf) if row.dd_drainaf is not None else None
    PSSI_depth = float(row.dd_draindi) if row.dd_draindi is not None else None
    measure = measure_func(SSI_distance, PSSI_distance)
    geometry = row.geometry_wkt

    postgres_insert_query = f"""
        INSERT INTO b_selection.b_{monitoringsjaar} (parcel_id, x, y, archetype, surface_level_ahn3, surface_level_ahn5, parcel_width, summer_stage, winter_stage, SSI_distance, SSI_depth,
        PSSI_distance, PSSI_depth, measure, geometry) VALUES (%(parcel_id)s,%(x)s,%(y)s,%(archetype)s,%(surface_level_ahn3)s,%(surface_level_ahn5)s,%(parcel_width)s,%(summer_stage)s,%(winter_stage)s,
                                            %(SSI_distance)s,%(SSI_depth)s,%(PSSI_distance)s,%(PSSI_depth)s,%(measure)s,%(geometry)s)
    """
    record_to_insert = {
        "parcel_id": str(parcel_id),
        "x": float(x),
        "y": float(y),
        "archetype": str(archetype),
        "surface_level_ahn3": float(surface_level_ahn3),
        "surface_level_ahn5": float(surface_level_ahn5),
        "parcel_width": float(parcel_width),
        "summer_stage": float(summer_stage),
        "winter_stage": float(winter_stage),
        "SSI_distance": SSI_distance,
        "SSI_depth": SSI_depth,
        "PSSI_distance": PSSI_distance,
        "PSSI_depth": PSSI_depth,
        "measure": measure,
        "geometry": geometry,
    }
    cursor.execute(postgres_insert_query, record_to_insert)
    connection.commit()
    count = cursor.rowcount

cursor.close()
connection.close()
print("PostgreSQL connection is closed")
