import geopandas as gpd
import pandas as pd
import numpy as np

POSSIBLE_ERRORS = {
    "is_insufficiently_covered": "Sum of Dekking_Veen and Dekking_Moerig is < 10%; ",
    "is_water": "Type is water; ",
    "is_buried_peat": "Buried peat soil; ",
    "is_petgat": "Archetype is a petgat, AVk; ",
    "is_invalid_surface_level": "Invalid surface level: should be in the range -8 - +3 (m NAP); ",
    "is_invalid_parcel_width": "Invalid parcel width. Should be between 5 and 300 meters; ",
    "is_invalid_summer_stage": "Summer stage should be below surface; ",
    "is_invalid_winter_stage": "Winter stage should be below surface; ",
    "is_invalid_summer_stage_freeboard": "Freeboard summer stage is < 2.m below surface level; ",
    "is_invalid_winter_stage_freeboard": "Freeboard winter stage is < 2.m below surface level; ",
    "is_invalid_PSSI_distance": "PSSI drain distance should be in range 3 - 10 (m); ",
    "is_invalid_SSI_distance": "SSI drain distance should be in range 3 - 10 (m); ",
    "is_invalid_SSI_depth": "Drain depth SSI: is defined as positive number below surface level and has to be within a depth of 1.2m; ",
    "is_invalid_PSSI_depth": "Drain depth PSSI: is defined as positive number below surface level and has to be within a depth of 1.2m; ",
}


def create_error_df(parcels_gdf: gpd.GeoDataFrame) -> pd.DataFrame:
    """
    Parameters
    ----------
    parcels_gdf : GeoDataFrame
        The GeoDataFrame containing parcel data.

    Returns
    -------
    errors : pandas.DataFrame
        A DataFrame with error flags for each parcel. The DataFrame includes the following columns:
            - Perceel_ID: Identifier of the parcel.
            - geometry: Geometry of the parcel.
            - is_insufficiently_covered: True if the sum of Dekking_Veen and Dekking_Moerig is less than 10.
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
    errors = pd.DataFrame(index=parcels_gdf.index)
    errors["is_insufficiently_covered"] = np.where(
        (parcels_gdf["Dekking_Veen"] + parcels_gdf["Dekking_Moerig"]) < 10, True, False
    )
    errors["is_water"] = np.where(parcels_gdf["type"] == "Water", True, False)
    errors["is_buried_peat"] = np.where(
        parcels_gdf["Archetype_Hoofd"].isin(["Rv01C", "AEm8", "Mv41C"]), True, False
    )
    errors["is_petgat"] = np.where(parcels_gdf["Archetype_Hoofd"] == "AVk", True, False)
    errors["is_invalid_surface_level"] = np.where(
        parcels_gdf["AHN3_Mediaan"].between(-8, 3), False, True
    )
    errors["is_invalid_parcel_width"] = np.where(
        parcels_gdf["Perceelsbreedte"].between(5, 300), False, True
    )
    errors["is_invalid_summer_stage"] = np.where(
        parcels_gdf["Peil_Gemiddeld_zp"] <= parcels_gdf["AHN3_Mediaan"], False, True
    )
    errors["is_invalid_winter_stage"] = np.where(
        (parcels_gdf["Peil_Gemiddeld_wp"] <= parcels_gdf["AHN3_Mediaan"]),
        False,
        True,
    )
    errors["is_invalid_summer_stage_freeboard"] = np.where(
        parcels_gdf["Peil_Gemiddeld_zp"] >= (parcels_gdf["AHN3_Mediaan"] - 2.0),
        False,
        True,
    )
    errors["is_invalid_winter_stage_freeboard"] = np.where(
        parcels_gdf["Peil_Gemiddeld_wp"] >= (parcels_gdf["AHN3_Mediaan"] - 2.0),
        False,
        True,
    )
    errors["is_invalid_PSSI_distance"] = np.where(
        (parcels_gdf["DD_Drainafstand"].isna())
        | (parcels_gdf["DD_Drainafstand"].between(3, 10)),
        False,
        True,
    )
    errors["is_invalid_SSI_distance"] = np.where(
        (parcels_gdf["OWD_Drainafstand"].isna())
        | (parcels_gdf["OWD_Drainafstand"].between(3, 10)),
        False,
        True,
    )
    errors["is_invalid_SSI_depth"] = np.where(
        (parcels_gdf["OWD_Draindiepte"].isna())
        | (parcels_gdf["OWD_Draindiepte"].between(0, 1.2)),
        False,
        True,
    )
    errors["is_invalid_PSSI_depth"] = np.where(
        (parcels_gdf["DD_Draindiepte"].isna())
        | (parcels_gdf["DD_Draindiepte"].between(0, 1.2)),
        False,
        True,
    )
    return errors


def filter_parcels(parcels_gdf: gpd.GeoDataFrame, errors: dict) -> gpd.GeoDataFrame:
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
    error_columns = [col for col in errors if col.startswith("is_")]
    invalid_mask = errors[error_columns].any(axis=1)
    valid_parcels = parcels_gdf.loc[~invalid_mask]

    invalid_parcels = parcels_gdf.loc[invalid_mask]["Perceel_ID"].copy().to_frame()
    invalid_parcels["error_messages"] = ""
    for error_type in POSSIBLE_ERRORS.keys():
        invalid_parcels.loc[errors[error_type], "error_messages"] += POSSIBLE_ERRORS[
            error_type
        ]

    return valid_parcels, invalid_parcels


def measure_func(SSI_dist, PSSI_dist):
    if not np.isnan(SSI_dist):
        measure = "SSI"
    elif not np.isnan(PSSI_dist):
        measure = "PSSI"
    else:
        measure = "ref"
    return measure


if __name__ == "__main__":
    # Load parcels data
    path_parcels = r"p:\11207812-somers-uitvoering\monitoring_2023\dataverzameling\updaten_AHN5.gpkg"
    parcels = gpd.read_file(path_parcels, layer="Percelen_2023_feb2025", crs=28992)

    # Add measure column (SSI, PSSI, or ref)
    parcels["measure"] = parcels.apply(
        lambda row: measure_func(row["OWD_Drainafstand"], row["DD_Drainafstand"]),
        axis=1,
    )

    # Get errors
    errors = create_error_df(parcels)
    valid_parcels, invalid_parcels = filter_parcels(parcels, errors)

    valid_parcels.to_file(
        r"P:/11207812-somers-uitvoering/monitoring_2023/shp_files/a_input/a_2023_testerik_fullset.shp"
    )
    invalid_parcels.to_csv(
        r"P:/11207812-somers-uitvoering/monitoring_2023/shp_files/a_input/error_lst_2023_testerik_fullset.csv"
    )
