import geopandas as gpd
import pandas as pd
import numpy as np

#%% Inputs
path_parcels = r'P:/11207812-somers-uitvoering/monitoring_2023/dataverzameling/updaten_AHN5.gpkg'

parcels_16_tot =  gpd.read_file(path_parcels, layer='Percelen_2016_feb2025', crs=28992)
parcels_22_tot =  gpd.read_file(path_parcels, layer='Percelen_2022_feb2025', crs=28992)
parcels_23_tot =  gpd.read_file(path_parcels, layer='Percelen_2023_feb2025', crs=28992)

#%% Min dekking 10%
parcels_16_tot = parcels_16_tot.loc[(parcels_16_tot['Dekking_Veen'] + (parcels_16_tot['Dekking_Moerig']) >= 10)]
parcels_22_tot = parcels_22_tot.loc[(parcels_22_tot['Dekking_Veen'] + (parcels_22_tot['Dekking_Moerig']) >= 10)]
parcels_23_tot = parcels_23_tot.loc[(parcels_23_tot['Dekking_Veen'] + (parcels_23_tot['Dekking_Moerig']) >= 10)]

#%% Remove water
def measure_func(SSI_dist, PSSI_dist):
    if not np.isnan(SSI_dist):
        measure = 'SSI'
    elif not np.isnan(PSSI_dist):
        measure = 'PSSI'
    else:
        measure = 'ref'
    return measure
        
        
error_lst = pd.DataFrame()

for n, row in parcels_22_tot.iterrows(): #check what record should be! And check order of variables in geopackage
    print(n)
    # Unpack row
    parcel_id = str(row.Perceel_ID)
    x = float(row.geometry.centroid.x)
    y = float(row.geometry.centroid.y)
    archetype = str(row.Archetype_Hoofd)
    bgt_type = str(row.type)
    surface_level_ahn3 = float(row.AHN3_Mediaan) if row.AHN3_Mediaan is not None else None
    surface_level_ahn5 = float(row.AHN5_Mediaan) if row.AHN5_Mediaan is not None else None
    parcel_width = float(row.Perceelsbreedte) if row.Perceelsbreedte is not None else None
    summer_stage= float(row.Peil_Gemiddeld_zp) if row.Peil_Gemiddeld_zp is not None else None
    winter_stage= float(row.Peil_Gemiddeld_wp) if row.Peil_Gemiddeld_wp is not None else None
    SSI_distance = float(row.OWD_Drainafstand) if row.OWD_Drainafstand is not None else None
    SSI_depth = float(row.OWD_Draindiepte) if row.OWD_Draindiepte is not None else None
    PSSI_distance = float(row.DD_Drainafstand) if row.DD_Drainafstand is not None else None
    PSSI_depth = float(row.DD_Draindiepte) if row.DD_Draindiepte is not None else None
    measure = measure_func(SSI_distance, PSSI_distance)
    geometry = row.geometry
    
    try:

        if any(
            value is None
            for value in [
                archetype, surface_level_ahn3, parcel_width, summer_stage, winter_stage
            ]
        ):
            error_lst.loc[parcel_id, 'error'] = "Missing a required parameter"
        if any(
            value is None or pd.isna(value)
            for value in [
                archetype, surface_level_ahn3, parcel_width, summer_stage, winter_stage
            ]
        ):
            error_lst.loc[parcel_id, 'error'] = "Missing a required parameter"
        elif archetype in ['Rv01C', 'AEm8', 'Mv41C']:
            error_lst.loc[parcel_id, 'error'] = "Burried peat soil"
        elif bgt_type == "Water":
            error_lst.loc[parcel_id, 'error'] = "Type is water"    
        elif archetype == 'AVk':
            error_lst.loc[parcel_id, 'error'] = "Archetype is a petgat, AVk"
        elif surface_level_ahn3 < -8 or surface_level_ahn3 > 3:
            error_lst.loc[parcel_id, 'error'] = "Invalid surface level: should be in the range -8 - +3 (m NAP)"
        elif parcel_width < 5 or parcel_width > 300:
            error_lst.loc[parcel_id, 'error'] = "Invalid parcel width. Should be between 5 and 300 meters"
        elif summer_stage > surface_level_ahn3:
            error_lst.loc[parcel_id, 'error'] = "Summer stage should be below surface"
        elif winter_stage > surface_level_ahn3:
            error_lst.loc[parcel_id, 'error'] = "Winter stage should be below surface"
        elif winter_stage < (surface_level_ahn3 - 2.0):
            error_lst.loc[parcel_id, 'error'] = "Winter stage is < 2.m below surface level"
        elif summer_stage < (surface_level_ahn3 - 2.0):
            error_lst.loc[parcel_id, 'error'] = "Summer stage is < 2.m below surface level"
        elif PSSI_distance is not None and (PSSI_distance > 10 or PSSI_distance < 3):
            error_lst.loc[parcel_id, 'error'] = "PSSI drain distance should be in range 3 - 10 (m)"
        elif SSI_distance is not None and (SSI_distance > 10 or SSI_distance < 3):
            error_lst.loc[parcel_id, 'error'] = "SSI drain distance should be in range 3 - 10 (m)"
        elif SSI_depth is not None and (SSI_depth < 0 or SSI_depth > 1.2):
            error_lst.loc[parcel_id, 'error'] = "Drain depth SSI: is defined as positive number below surface level and has to be within a depth of 1.2m"
        elif PSSI_depth is not None and (PSSI_depth < 0 or PSSI_depth > 1.2):
            error_lst.loc[parcel_id, 'error'] = "Drain depth PSSI: is defined as positive number below surface level and has to be within a depth of 1.2m"    
    except Exception as e:
        error_lst.loc[parcel_id] = f"Error processing row: {e}"

# Keep rows without errors
valid_rows = parcels_22_tot.loc[~parcels_22_tot['Perceel_ID'].isin(error_lst.index)]

#%% Saving
valid_rows.to_file(r'P:/11207812-somers-uitvoering/monitoring_2023/shp_files/a_input/a_2022.shp')
error_lst.to_csv(r'P:/11207812-somers-uitvoering/monitoring_2023/shp_files/a_input/error_lst_2022.csv')