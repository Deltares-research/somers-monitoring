import pandas as pd
import warnings
import configparser
import psycopg2
import geopandas as gpd
import sys
from shapely.geometry import Point
import os
import numpy as np
os.environ['USE_PYGEOS'] = '0'
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
warnings.filterwarnings("ignore", category=DeprecationWarning)

import somers
from somers import to_linux
sys.path.append(to_linux(r'P:/11207812-somers-uitvoering/monitoring_2023/src/'))


#%%

b_monitoringsjaar = gpd.read_file(to_linux(r'C:/Users/nougu_la/Desktop/test_percelen.gpkg'))
b_monitoringsjaar = b_monitoringsjaar.set_index('PARCEL_ID')

parcels = b_monitoringsjaar.copy()

#%%
path_names = pd.read_csv(to_linux(r"P:/11207812-somers-ontwikkeling/somers_v2.1_input/data/pathnames_monitoring.csv"), index_col = 0 ,delimiter = ";").squeeze().to_dict()
path_names  = {key: to_linux(value) for key, value in path_names.items()}

def find_weerjaar(x, y):
    point = Point(x, y)
    weerregios = gpd.read_file(to_linux(r'P:/11207812-somers-ontwikkeling/somers_v2.1_input/data/meteo/weather_regions.shp'))
    weerregio = weerregios[weerregios.contains(point)]['weather_rg'].item()
    if weerregio == 'southwest':
        return '2015' #2015 noramly
    elif weerregio == 'northeast':
        return '2013' #2013 normaly
    
        
def find_ditch_depth(surface_level: float, ditch_stage: float, default_ditch_depth: float=1.2, min_ditch_depth: float=0.4):
    default_water_depth = default_ditch_depth - (surface_level - ditch_stage)
    if default_water_depth < min_ditch_depth:
        ditch_depth = surface_level - (ditch_stage - min_ditch_depth)
    else:
        ditch_depth = default_ditch_depth
    return ditch_depth   
    

def set_max_freeboard(surface_level: float, ditch_stage: float, max_freeboard: float=1.5):
    if (surface_level - ditch_stage) > max_freeboard:
        return (surface_level - max_freeboard)
    else: 
        return ditch_stage

def set_min_freeboard(surface_level: float, ditch_stage: float, min_freeboard: float=0.051):
    if (surface_level - ditch_stage) < min_freeboard:
        return (surface_level - min_freeboard)
    else: 
        return ditch_stage
    
def adjust_freeboard(surface_level: float, ditch_stage: float):
    ditch_stage_1 = set_max_freeboard(surface_level, ditch_stage)
    ditch_stage_2 = set_min_freeboard(surface_level, ditch_stage_1)
    
    return ditch_stage_2


def run_somers_ref(parcel_id, surface_level, parcel_width, winter_stage, summer_stage, x, y, archetype):
    parameter_combs = pd.read_csv(to_linux(r'P:/11207812-somers-ontwikkeling/somers_v2.1_development/modflow_calibration/4-output/somers_2_1_parameter_selection.csv'), index_col = 0)
    
    try:
        output = somers.run_somers(
            ## model settings
            name = str(parcel_id),
            path_names = path_names,
            input_dir=to_linux(r'P:/11207812-somers-uitvoering/monitoring_2023/run_dir/3-input'),
            output_dir=to_linux(r'P:/11207812-somers-uitvoering/monitoring_2023/run_dir/4-output'),
            working_dir=to_linux(r'P:/11207812-somers-uitvoering/monitoring_2023/run_dir/2-interim'),
            #modflow_executable=to_linux(r'P:/11207812-somers-uitvoering/software/mf6_compiled_rm_v3'), #linux
            modflow_executable=to_linux(r'P:/11207812-somers-uitvoering/software/mf6.exe'), #local test
            dimension = '2D',
            method_groundwater = 'modflow', 
            method_aquifer = 'flux',
            method_soil_moisture = 'dynamic',
            method_soil_temperature = 'fft',
            method_aap_wfps = 'somers', 
            method_aap_temp = 'somers', 
            method_br = 'somers_v2',
            parameter_file = parameter_combs,
        
            ## parcel properties
            parcel_type = 'ref',
            x=x,
            y=y,
            surface_level=float(surface_level),
            parcel_width=float(parcel_width),
            start_date=f'01-01-{find_weerjaar(x, y)}',
            end_date=f'31-12-{find_weerjaar(x, y)}',
        
            ## ditch parameters
            ditch_depth=find_ditch_depth(surface_level, np.min([winter_stage, summer_stage])),
            winter_stage = float(winter_stage),
            summer_stage = float(summer_stage),
            
            ## preset values
            preset_soilcode = archetype, 
            
            ## output
            save_input = False,
            save_modflow_runs = False,
            save_phreatic_head = False,
            save_soil_moisture = False,
            save_soil_temperature = False,
            save_emission= False,
            
            )
                
        output_median = output['emission_be'].sum(['time','depth','x']).median('runs').item()

        
        return (output_median)

        
    except (psycopg2.errors.UniqueViolation, ValueError, FloatingPointError, KeyError, FileNotFoundError, IndexError, TypeError, np.core._exceptions._ArrayMemoryError, PermissionError) as error:
        
        with open(to_linux(r"P:/11207812-somers-uitvoering/monitoring_2023/shp_files/error_lst_ref.txt"), "a") as file:
            file.write('parcel: '+str(parcel_id)+ "; " + str(error) + "\n")
            
        return (None)
#%%
print(f'running n={len(parcels)} parcels')

parcel_id_lst = []
output_median_lst = []

runs = pd.DataFrame()    
for parcel_id, parcel in parcels.iterrows():   
    print(parcel_id)
    
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
        

    print('Run ref model sim')
    #continue
    output_median = run_somers_ref(
        parcel_id=parcel_id,
        x=x,
        y=y,
        archetype=archetype,
        parcel_width=parcel_width,
        surface_level=surface_level,
        summer_stage=summer_stage,
        winter_stage=winter_stage)

    parcel_id_lst.append(parcel_id)
    output_median_lst.append(output_median)
    
normal_run = pd.DataFrame()
normal_run['parcel_id'] = parcel_id_lst
normal_run = normal_run.set_index(normal_run.parcel_id, drop = True)
normal_run['emission'] = output_median_lst   
normal_run['parcel_width'] = parcels.PARCEL_WID
normal_run['parcel_area'] = parcels.geometry.area #m2
normal_run['parcel_length'] = normal_run['parcel_area'] / normal_run['parcel_width']

normal_run.to_csv(r'C:/Users/nougu_la/Desktop/tnormal_run.csv')
#%% bigger width

parcel_id_wide_lst = []
output_median_wide_lst = []
parcel_width_lst = []

for parcel_id, parcel in parcels.iterrows():   
    print(parcel_id)
    
    ## load input
    x = parcel.X
    y = parcel.Y
    archetype = parcel.ARCHETYPE
    parcel_width = parcel.PARCEL_WID + 7
    surface_level = parcel.SURFACE_LE #AHN3
    summer_stage = adjust_freeboard(surface_level, parcel.SUMMER_STA)
    winter_stage = adjust_freeboard(surface_level, parcel.WINTER_STA)
    ssi_distance = parcel.SSI_DISTAN
    ssi_depth = parcel.SSI_DEPTH
    pssi_distance = parcel.PSSI_DISTA
    pssi_depth = parcel.PSSI_DEPTH
    measure = parcel.MEASURE
    geometry = parcel.geometry
        

    print('Run ref model sim')
    #continue
    output_median = run_somers_ref(
        parcel_id=parcel_id,
        x=x,
        y=y,
        archetype=archetype,
        parcel_width=parcel_width,
        surface_level=surface_level,
        summer_stage=summer_stage,
        winter_stage=winter_stage)

    parcel_id_wide_lst.append(parcel_id)
    output_median_wide_lst.append(output_median)
    parcel_width_lst.append(parcel_width)
    
wide_run = pd.DataFrame()
wide_run['parcel_id'] = parcel_id_wide_lst
wide_run = wide_run.set_index(wide_run.parcel_id, drop = True)
wide_run['emission'] = output_median_wide_lst   
wide_run['parcel_width'] = parcel_width_lst
wide_run['parcel_length'] = normal_run['parcel_length']
wide_run['parcel_area'] = wide_run['parcel_width']  * wide_run['parcel_length']

wide_run.to_csv(r'C:/Users/nougu_la/Desktop/wide_run.csv')
