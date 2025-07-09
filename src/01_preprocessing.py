import geopandas as gpd
import pandas as pd
import numpy as np

# %% Inputs
path_parcels = (
    r"p:\11207812-somers-uitvoering\monitoring_2023\dataverzameling\test_percelen.gpkg"
)
path_parcels = r"p:\11207812-somers-uitvoering\monitoring_2023\shp_files\a_input\a_2023_testerik.shp"
parcels_16_tot = gpd.read_file(path_parcels)

parcels_16_tot = gpd.read_file(
    path_parcels, layer="updaten_ahn5__percelen_2016_feb2025", crs=28992
)
# parcels_22_tot =  gpd.read_file(path_parcels, layer='Percelen_2022_feb2025', crs=28992)
# parcels_23_tot =  gpd.read_file(path_parcels, layer='Percelen_2023_feb2025', crs=28992)

# %% Min dekking 10%
parcels_16_tot = parcels_16_tot.loc[
    (parcels_16_tot["Dekking_Veen"] + (parcels_16_tot["Dekking_Moerig"]) >= 10)
]
# parcels_22_tot = parcels_22_tot.loc[(parcels_22_tot['Dekking_Veen'] + (parcels_22_tot['Dekking_Moerig']) >= 10)]
# parcels_23_tot = parcels_23_tot.loc[(parcels_23_tot['Dekking_Veen'] + (parcels_23_tot['Dekking_Moerig']) >= 10)]

ERROR_MESSAGES = {
    "is_water": "Type is water",
    "is_buried_peat": "Buried peat soil",
    "is_petgat": "Archetype is a petgat, AVk",
    "is_invalid_surface_level": "Invalid surface level: should be in the range -8 - +3 (m NAP)",
    "is_invalid_parcel_width": "Invalid parcel width. Should be between 5 and 300 meters",
    "is_invalid_summer_stage": "Summer stage should be below surface",
    "is_invalid_winter_stage": "Winter stage should be below surface",
    "is_invalid_summer_stage_freeboard": "Freeboard summer stage is < 2.m below surface level",
    "is_invalid_winter_stage_freeboard": "Freeboard winter stage is < 2.m below surface level",
    "is_invalid_PSSI_distance": "PSSI drain distance should be in range 3 - 10 (m)",
    "is_invalid_SSI_distance": "SSI drain distance should be in range 3 - 10 (m)",
    "is_invalid_SSI_depth": "Drain depth SSI: is defined as positive number below surface level and has to be within a depth of 1.2m",
    "is_invalid_PSSI_depth": "Drain depth PSSI: is defined as positive number below surface level and has to be within a depth of 1.2m",
}


def create_error_df(parcels_gdf):
    """
    Parameters
    ----------
    parcels_gdf : GeoDataFrame
        The GeoDataFrame containing parcel data.
    Returns
    -------
    error_df : DataFrame
        A DataFrame with error flags for each parcel. The DataFrame includes the following columns:
            - Perceel_ID: Identifier of the parcel.
            - geometry: Geometry of the parcel.
            - is_water: True if the parcel type is 'Water'.
            - is_buried_peat: True if the parcel archetype is one of ['Rv01C', 'AEm8', 'Mv41C'].
            - is_petgat: True if the parcel archetype is 'AVk'.
            - is_invalid_surface_level: True if AHN3_Mediaan is not between -8 and 3.
            - is_invalid_parcel_width: True if Perceelsbreedte is not between 5 and 300.
            - is_invalid_summer_stage: True if Peil_Gemiddeld_zp is greater than AHN3_Mediaan.
            - is_invalid_winter_stage: True if Peil_Gemiddeld_wp is not between (AHN3_Mediaan - 2.0) and AHN3_Mediaan.
            - is_invalid_PSSI_distance: True if DD_Drainafstand is not between 3 and 10 (when not null).
            - is_invalid_SSI_distance: True if OWD_Drainafstand is not between 3 and 10 (when not null).
            - is_invalid_SSI_depth: True if OWD_Draindiepte is not between 0 and 1.2 (when not null).
            - is_invalid_PSSI_depth: True if DD_Draindiepte is not between 0 and 1.2 (when not null).
    Notes
    -----
    Each error flag column is a boolean indicating whether the corresponding condition is met for each parcel.
    """
    error_df = parcels_gdf[["Perceel_ID", "geometry"]].copy()
    error_df["is_water"] = np.where(parcels_gdf["type"] == "Water", True, False)
    error_df["is_buried_peat"] = np.where(
        parcels_gdf["Archetype_Hoofd"].isin(["Rv01C", "AEm8", "Mv41C"]), True, False
    )
    error_df["is_petgat"] = np.where(
        parcels_gdf["Archetype_Hoofd"] == "AVk", True, False
    )
    error_df["is_invalid_surface_level"] = np.where(
        parcels_gdf["AHN3_Mediaan"].between(-8, 3), False, True
    )
    error_df["is_invalid_parcel_width"] = np.where(
        parcels_gdf["Perceelsbreedte"].between(5, 300), False, True
    )
    error_df["is_invalid_summer_stage"] = np.where(
        parcels_gdf["Peil_Gemiddeld_zp"] <= parcels_gdf["AHN3_Mediaan"], False, True
    )
    error_df["is_invalid_winter_stage"] = np.where(
        (parcels_gdf["Peil_Gemiddeld_wp"] <= parcels_gdf["AHN3_Mediaan"]),
        False,
        True,
    )
    error_df["is_invalid_summer_stage_freeboard"] = np.where(
        parcels_gdf["Peil_Gemiddeld_zp"] >= (parcels_gdf["AHN3_Mediaan"] - 2.0),
        False,
        True,
    )
    error_df["is_invalid_winter_stage_freeboard"] = np.where(
        parcels_gdf["Peil_Gemiddeld_wp"] >= (parcels_gdf["AHN3_Mediaan"] - 2.0),
        False,
        True,
    )
    error_df["is_invalid_PSSI_distance"] = np.where(
        (parcels_gdf["DD_Drainafstand"].isna())
        | (parcels_gdf["DD_Drainafstand"].between(3, 10)),
        False,
        True,
    )
    error_df["is_invalid_SSI_distance"] = np.where(
        (parcels_gdf["OWD_Drainafstand"].isna())
        | (parcels_gdf["OWD_Drainafstand"].between(3, 10)),
        False,
        True,
    )
    error_df["is_invalid_SSI_depth"] = np.where(
        (parcels_gdf["OWD_Draindiepte"].isna())
        | (parcels_gdf["OWD_Draindiepte"].between(0, 1.2)),
        False,
        True,
    )
    error_df["is_invalid_PSSI_depth"] = np.where(
        (parcels_gdf["DD_Draindiepte"].isna())
        | (parcels_gdf["DD_Draindiepte"].between(0, 1.2)),
        False,
        True,
    )
    return error_df


def filter_parcels(parcels_gdf, error_df):
    """
    Parameters
    ----------
    parcels_gdf : GeoDataFrame
        The GeoDataFrame containing parcel data.
    error_df : DataFrame
        The DataFrame with error flags for each parcel.

    Returns
    -------
    valid_parcels : GeoDataFrame
        A GeoDataFrame containing only the parcels that do not have any errors.
    """
    # Combine all error columns (those starting with 'is_') to filter out any parcels with errors
    error_columns = [col for col in error_df.columns if col.startswith("is_")]
    valid_mask = ~error_df[error_columns].any(axis=1)
    valid_parcels = parcels_gdf.loc[valid_mask]
    return valid_parcels


def measure_func(SSI_dist, PSSI_dist):
    if not np.isnan(SSI_dist):
        measure = "SSI"
    elif not np.isnan(PSSI_dist):
        measure = "PSSI"
    else:
        measure = "ref"
    return measure


parcels_16_tot["measure"] = parcels_16_tot.apply(
    lambda row: measure_func(row["OWD_Drainafstand"], row["DD_Drainafstand"]),
    axis=1,
)
error_df = create_error_df(parcels_16_tot)
valid_parcels = filter_parcels(parcels_16_tot, error_df)

## TODO: error logging

# %% BELOW IS OLD CODE, NOT USED ANYMORE

# %% Remove water


error_lst = pd.DataFrame()

for (
    n,
    row,
) in (
    parcels_16_tot.iterrows()
):  # check what record should be! And check order of variables in geopackage
    print(n)
    # Unpack row
    parcel_id = str(row.Perceel_ID)
    x = float(row.geometry.centroid.x)
    y = float(row.geometry.centroid.y)
    archetype = str(row.Archetype_Hoofd)
    bgt_type = str(row.type)
    surface_level_ahn3 = (
        float(row.AHN3_Mediaan) if row.AHN3_Mediaan is not None else None
    )
    surface_level_ahn5 = (
        float(row.AHN5_Mediaan) if row.AHN5_Mediaan is not None else None
    )
    parcel_width = (
        float(row.Perceelsbreedte) if row.Perceelsbreedte is not None else None
    )
    summer_stage = (
        float(row.Peil_Gemiddeld_zp) if row.Peil_Gemiddeld_zp is not None else None
    )
    winter_stage = (
        float(row.Peil_Gemiddeld_wp) if row.Peil_Gemiddeld_wp is not None else None
    )
    SSI_distance = (
        float(row.OWD_Drainafstand) if row.OWD_Drainafstand is not None else None
    )
    SSI_depth = float(row.OWD_Draindiepte) if row.OWD_Draindiepte is not None else None
    PSSI_distance = (
        float(row.DD_Drainafstand) if row.DD_Drainafstand is not None else None
    )
    PSSI_depth = float(row.DD_Draindiepte) if row.DD_Draindiepte is not None else None
    measure = measure_func(SSI_distance, PSSI_distance)
    geometry = row.geometry

    try:

        if any(
            value is None or pd.isna(value)
            for value in [
                archetype,
                surface_level_ahn3,
                parcel_width,
                summer_stage,
                winter_stage,
            ]
        ):
            error_lst.loc[parcel_id, "error"] = "Missing a required parameter"
        elif archetype in ["Rv01C", "AEm8", "Mv41C"]:
            error_lst.loc[parcel_id, "error"] = "Burried peat soil"
        elif bgt_type == "Water":
            error_lst.loc[parcel_id, "error"] = "Type is water"
        elif archetype == "AVk":
            error_lst.loc[parcel_id, "error"] = "Archetype is a petgat, AVk"
        elif surface_level_ahn3 < -8 or surface_level_ahn3 > 3:
            error_lst.loc[parcel_id, "error"] = (
                "Invalid surface level: should be in the range -8 - +3 (m NAP)"
            )
        elif parcel_width < 5 or parcel_width > 300:
            error_lst.loc[parcel_id, "error"] = (
                "Invalid parcel width. Should be between 5 and 300 meters"
            )
        elif summer_stage > surface_level_ahn3:
            error_lst.loc[parcel_id, "error"] = "Summer stage should be below surface"
        elif winter_stage > surface_level_ahn3:
            error_lst.loc[parcel_id, "error"] = "Winter stage should be below surface"
        elif winter_stage < (surface_level_ahn3 - 2.0):
            error_lst.loc[parcel_id, "error"] = (
                "Winter stage is < 2.m below surface level"
            )
        elif summer_stage < (surface_level_ahn3 - 2.0):
            error_lst.loc[parcel_id, "error"] = (
                "Summer stage is < 2.m below surface level"
            )
        elif PSSI_distance is not None and (PSSI_distance > 10 or PSSI_distance < 3):
            error_lst.loc[parcel_id, "error"] = (
                "PSSI drain distance should be in range 3 - 10 (m)"
            )
        elif SSI_distance is not None and (SSI_distance > 10 or SSI_distance < 3):
            error_lst.loc[parcel_id, "error"] = (
                "SSI drain distance should be in range 3 - 10 (m)"
            )
        elif SSI_depth is not None and (SSI_depth < 0 or SSI_depth > 1.2):
            error_lst.loc[parcel_id, "error"] = (
                "Drain depth SSI: is defined as positive number below surface level and has to be within a depth of 1.2m"
            )
        elif PSSI_depth is not None and (PSSI_depth < 0 or PSSI_depth > 1.2):
            error_lst.loc[parcel_id, "error"] = (
                "Drain depth PSSI: is defined as positive number below surface level and has to be within a depth of 1.2m"
            )
    except Exception as e:
        error_lst.loc[parcel_id] = f"Error processing row: {e}"

# Keep rows without errors
valid_rows = parcels_16_tot.loc[~parcels_16_tot["Perceel_ID"].isin(error_lst.index)]

# %% Saving
valid_rows.to_file(
    r"P:/11207812-somers-uitvoering/monitoring_2023/shp_files/a_input/a_2023_testerik.shp"
)
error_lst.to_csv(
    r"P:/11207812-somers-uitvoering/monitoring_2023/shp_files/a_input/error_lst_2023_testerik.csv"
)
