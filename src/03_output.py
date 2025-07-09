import pandas as pd
import warnings
import configparser
import psycopg2
import sys
import os
import numpy as np

os.environ["USE_PYGEOS"] = "0"
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
warnings.filterwarnings("ignore", category=DeprecationWarning)

from somers import to_linux

sys.path.append(to_linux(r"P:/11207812-somers-uitvoering/monitoring_2023/src/"))
from output_functions import (
    adjust_freeboard,
    run_somers_ref,
    run_somers_ssi,
    run_somers_pssi,
)

# %% Database connection
print("Starting initial connection to PostgreSQL")

# Read the variables from the command-line arguments
variables = sys.argv[1:]
# variables = [2023, 0, 2000] ## DEZE MOET JE COMMENTEN NA HET PROBEREN
# print(sys.argv[1:])
monitorings_jaar_idx = int(variables[0])
start_idx = int(variables[1])
end_idx = int(variables[2])

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

monitoringsjaar = monitorings_jaar_idx
# %%

sql_create_table = f"""CREATE TABLE IF NOT EXISTS c_output.c_{monitoringsjaar}(
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
    succesfull_runs INT NOT NULL,
    measure TEXT,
    geometry GEOMETRY);

    ALTER TABLE c_output.c_{monitoringsjaar} OWNER TO somers"""

cursor.execute(sql_create_table)
connection.commit()

# Download B scheme
sql_select_query = f"""SELECT parcel_id, x, y, archetype, surface_level_ahn3, surface_level_ahn5, parcel_width, 
    summer_stage, winter_stage, ssi_distance, ssi_depth, pssi_distance, pssi_depth, measure, geometry 
    FROM b_selection.b_{monitoringsjaar} ORDER BY parcel_id"""

cursor.execute(sql_select_query)
record = cursor.fetchall()

connection.commit()

columns = {
    "parcel_id": str,
    "x": float,
    "y": float,
    "archetype": str,
    "surface_level_ahn3": float,
    "surface_level_ahn5": float,
    "parcel_width": float,
    "summer_stage": float,
    "winter_stage": float,
    "ssi_distance": float,
    "ssi_depth": float,
    "pssi_distance": float,
    "pssi_depth": float,
    "measure": str,
    "geometry": object,
}

b_monitoringsjaar = (
    pd.DataFrame(record, columns=columns.keys()).astype(columns).set_index("parcel_id")
)

# %%Check which parcels have already been calculated
sql_completed_parcels = f"""SELECT DISTINCT parcel_id 
       FROM c_output.c_{monitoringsjaar}"""
cursor.execute(sql_completed_parcels)
completed_parcels = [str(i[0]) for i in cursor.fetchall()]

cursor.close()
connection.close()

ssi_depth_mean = b_monitoringsjaar.ssi_depth.mean()

print("Initial PostgreSQL connection is closed")
# %%
print("Starting model run connection to PostgreSQL")

skip_lst = [
    "RP-L-3077305",
    "RP-L-3089835",
    "RP-L-3296055",
    "RP-L-3417857",
    "RP-L-3305432",
    "RP-L-3417857",
    "RP-L-3417726",
    "RP-L-3073591",
    "RP-L-3353583",
    "RP-L-3453667",
    "RP-L-3129115",
    "RP-H-2297354",
    "RP-L-3305970 ",
]

parcels = b_monitoringsjaar[~b_monitoringsjaar.index.isin(completed_parcels)]
parcels = parcels[~parcels.index.isin(skip_lst)]
parcels = parcels.iloc[start_idx : end_idx + 1]
# %%
print(f"running n={len(parcels)} parcels")


runs = pd.DataFrame()
for parcel_id, parcel in parcels.iterrows():
    print(parcel_id)

    # make connection
    connection = psycopg2.connect(
        host=host, database=database, user=user, password=password, port=port
    )
    cursor = connection.cursor()

    ## load input
    x = parcel.x
    y = parcel.y
    archetype = parcel.archetype
    parcel_width = parcel.parcel_width
    surface_level = parcel.surface_level_ahn5
    summer_stage = adjust_freeboard(surface_level, parcel.summer_stage)
    winter_stage = adjust_freeboard(surface_level, parcel.winter_stage)
    ssi_distance = parcel.ssi_distance
    ssi_depth = parcel.ssi_depth
    pssi_distance = parcel.pssi_distance
    pssi_depth = parcel.pssi_depth
    measure = parcel.measure
    geometry = parcel.geometry

    if parcel.measure == "SSI":
        print("Run SSI model sim")
        # continue

        if np.isnan(parcel.ssi_depth):
            ssi_depth = ssi_depth_mean
        else:
            ssi_depth = ssi_depth

        (
            output_min,
            output_median,
            output_max,
            output_10_perc,
            output_90_perc,
            min_summer_gwl,
            min_summer3_gwl,
            min_winter_gwl,
            min_year_gwl,
            min_rlg_gwl,
            min_rhg_gwl,
            mean_summer_gwl,
            mean_summer3_gwl,
            mean_winter_gwl,
            mean_year_gwl,
            mean_rlg_gwl,
            mean_rhg_gwl,
            max_summer_gwl,
            max_summer3_gwl,
            max_winter_gwl,
            max_year_gwl,
            max_rlg_gwl,
            max_rhg_gwl,
            nb_runs,
        ) = run_somers_ssi(
            parcel_id=parcel_id,
            x=x,
            y=y,
            archetype=archetype,
            parcel_width=parcel_width,
            surface_level=surface_level,
            summer_stage=summer_stage,
            winter_stage=winter_stage,
            ssi_distance=ssi_distance,
            ssi_depth=ssi_depth,
        )

    elif parcel.measure == "PSSI":
        print("Run PSSI model sim")
        # continue
        (
            output_min,
            output_median,
            output_max,
            output_10_perc,
            output_90_perc,
            min_summer_gwl,
            min_summer3_gwl,
            min_winter_gwl,
            min_year_gwl,
            min_rlg_gwl,
            min_rhg_gwl,
            mean_summer_gwl,
            mean_summer3_gwl,
            mean_winter_gwl,
            mean_year_gwl,
            mean_rlg_gwl,
            mean_rhg_gwl,
            max_summer_gwl,
            max_summer3_gwl,
            max_winter_gwl,
            max_year_gwl,
            max_rlg_gwl,
            max_rhg_gwl,
            nb_runs,
        ) = run_somers_pssi(
            parcel_id=parcel_id,
            x=x,
            y=y,
            archetype=archetype,
            parcel_width=parcel_width,
            surface_level=surface_level,
            summer_stage=summer_stage,
            winter_stage=winter_stage,
            pssi_distance=pssi_distance,
            pssi_depth=pssi_depth,
        )

    elif parcel.measure == "ref":
        print("Run ref model sim")
        # continue
        (
            output_min,
            output_median,
            output_max,
            output_10_perc,
            output_90_perc,
            min_summer_gwl,
            min_summer3_gwl,
            min_winter_gwl,
            min_year_gwl,
            min_rlg_gwl,
            min_rhg_gwl,
            mean_summer_gwl,
            mean_summer3_gwl,
            mean_winter_gwl,
            mean_year_gwl,
            mean_rlg_gwl,
            mean_rhg_gwl,
            max_summer_gwl,
            max_summer3_gwl,
            max_winter_gwl,
            max_year_gwl,
            max_rlg_gwl,
            max_rhg_gwl,
            nb_runs,
        ) = run_somers_ref(
            parcel_id=parcel_id,
            x=x,
            y=y,
            archetype=archetype,
            parcel_width=parcel_width,
            surface_level=surface_level,
            summer_stage=summer_stage,
            winter_stage=winter_stage,
        )

    if output_min is None:
        continue

    postgres_insert_query = f"""INSERT INTO c_output.c_{monitoringsjaar} (parcel_id, minimum, median, maximum, perc_10, perc_90, 
        min_summer_gwl,min_summer3_gwl,min_winter_gwl,min_year_gwl,min_rlg_gwl,min_rhg_gwl,
        mean_summer_gwl,mean_summer3_gwl,mean_winter_gwl,mean_year_gwl,mean_rlg_gwl,mean_rhg_gwl,
        max_summer_gwl,max_summer3_gwl,max_winter_gwl,max_year_gwl,max_rlg_gwl,max_rhg_gwl,
        succesfull_runs, measure, geometry)
        VALUES (%(parcel_id)s, %(minimum)s, %(median)s, %(maximum)s, %(perc_10)s, %(perc_90)s, 
        %(min_summer_gwl)s, %(min_summer3_gwl)s, %(min_winter_gwl)s, %(min_year_gwl)s, %(min_rlg_gwl)s, %(min_rhg_gwl)s, 
        %(mean_summer_gwl)s, %(mean_summer3_gwl)s, %(mean_winter_gwl)s, %(mean_year_gwl)s, %(mean_rlg_gwl)s, %(mean_rhg_gwl)s, 
        %(max_summer_gwl)s, %(max_summer3_gwl)s, %(max_winter_gwl)s, %(max_year_gwl)s, %(max_rlg_gwl)s, %(max_rhg_gwl)s,
        %(succesfull_runs)s, %(measure)s, %(geometry)s)
        ON CONFLICT (parcel_id) DO UPDATE 
        SET minimum = EXCLUDED.minimum, median = EXCLUDED.median, maximum = EXCLUDED.maximum, perc_10 = EXCLUDED.perc_10, perc_90 = EXCLUDED.perc_90, 
        min_summer_gwl = Excluded.min_summer_gwl, min_summer3_gwl = Excluded.min_summer3_gwl, min_winter_gwl = Excluded.min_winter_gwl, min_year_gwl = Excluded.min_year_gwl, min_rlg_gwl = Excluded.min_rlg_gwl, min_rhg_gwl = Excluded.min_rhg_gwl,
        mean_summer_gwl = Excluded.mean_summer_gwl, mean_summer3_gwl = Excluded.mean_summer3_gwl, mean_winter_gwl = Excluded.mean_winter_gwl, mean_year_gwl = Excluded.mean_year_gwl, mean_rlg_gwl = Excluded.mean_rlg_gwl, mean_rhg_gwl = Excluded.mean_rhg_gwl,
        max_summer_gwl = Excluded.max_summer_gwl, max_summer3_gwl = Excluded.max_summer3_gwl, max_winter_gwl = Excluded.max_winter_gwl, max_year_gwl = Excluded.max_year_gwl, max_rlg_gwl = Excluded.max_rlg_gwl, max_rhg_gwl = Excluded.max_rhg_gwl,
        succesfull_runs = EXCLUDED.succesfull_runs, measure = EXCLUDED.measure, geometry = EXCLUDED.geometry"""
    record_to_insert = {
        "parcel_id": str(parcel_id),
        "minimum": float(output_min),
        "median": float(output_median),
        "maximum": float(output_max),
        "perc_10": float(output_10_perc),
        "perc_90": float(output_90_perc),
        "min_summer_gwl": float(min_summer_gwl),
        "min_summer3_gwl": float(min_summer3_gwl),
        "min_winter_gwl": float(min_winter_gwl),
        "min_year_gwl": float(min_year_gwl),
        "min_rlg_gwl": float(min_rlg_gwl),
        "min_rhg_gwl": float(min_rhg_gwl),
        "mean_summer_gwl": float(mean_summer_gwl),
        "mean_summer3_gwl": float(mean_summer3_gwl),
        "mean_winter_gwl": float(mean_winter_gwl),
        "mean_year_gwl": float(mean_year_gwl),
        "mean_rlg_gwl": float(mean_rlg_gwl),
        "mean_rhg_gwl": float(mean_rhg_gwl),
        "max_summer_gwl": float(max_summer_gwl),
        "max_summer3_gwl": float(max_summer3_gwl),
        "max_winter_gwl": float(max_winter_gwl),
        "max_year_gwl": float(max_year_gwl),
        "max_rlg_gwl": float(max_rlg_gwl),
        "max_rhg_gwl": float(max_rhg_gwl),
        "succesfull_runs": int(nb_runs),
        "measure": str(measure),
        "geometry": geometry,
    }
    cursor.execute(postgres_insert_query, record_to_insert)

    # output is in kg co2/ha/jaar
    connection.commit()
    count = cursor.rowcount
    print(count, "Record inserted successfully into mobile table")

    cursor.close()
    connection.close()
    print("PostgreSQL connection for single run is closed")

print("PostgreSQL connection for all runs is closed")
