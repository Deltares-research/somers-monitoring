import numpy as np
import pandas as pd
import psycopg2
import xarray as xr
import geopandas as gpd
from shapely.geometry import Point
from datetime import timedelta

import os
os.environ['USE_PYGEOS'] = '0'
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

import somers
from somers import to_linux

#%% SOMERS functions

path_names = pd.read_csv(to_linux(r"P:/11207812-somers-ontwikkeling/somers_v2.1_input/data/pathnames_monitoring.csv"), index_col = 0 ,delimiter = ";").squeeze().to_dict()
path_names  = {key: to_linux(value) for key, value in path_names.items()}

#functie scrhijven voor parametercombs bestand afhankelijk van archetype
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

    
def find_ditch_depth(surface_level: float, ditch_stage: float, default_ditch_depth: float=1.2, min_ditch_depth: float=0.4):
    default_water_depth = default_ditch_depth - (surface_level - ditch_stage)
    if default_water_depth < min_ditch_depth:
        ditch_depth = surface_level - (ditch_stage - min_ditch_depth)
    else:
        ditch_depth = default_ditch_depth
    return ditch_depth   
    
def find_weerjaar(x, y):
    point = Point(x, y)
    weerregios = gpd.read_file(to_linux(r'P:/11207812-somers-ontwikkeling/somers_v2.1_input/data/meteo/weather_regions.shp'))
    weerregio = weerregios[weerregios.contains(point)]['weather_rg'].item()
    if weerregio == 'southwest':
        return '2015'
    elif weerregio == 'northeast':
        return '2013'

def load_recharge(
    x: float,
    y: float,
    path_names: dict,
):
    weather_year = find_weerjaar(x, y)
    start_date = pd.to_datetime(f'01-01-{weather_year}', format = '%d-%m-%Y')
    end_date = pd.to_datetime(f'31-12-{weather_year}', format = '%d-%m-%Y')
    
    lhm_recharge = xr.open_dataarray(path_names["path_lhm_recharge"])
    recharge_at_xy = lhm_recharge.sel(x=x, y=y, method="nearest")
    recharge_series = recharge_at_xy.sel(
        time=slice(start_date, end_date)
            )/ 1000 
    return recharge_series.drop_vars(['x','y','dx','dy']).to_dataframe('recharge (m/d)')

def set_pssi_scenario(
        scenario: str,
        x: float, 
        y: float,
        surface_level: float,
        summer_stage: float,
        path_names: dict,
        infiltration_threshold_high: float=0.000,
        infiltration_threshold_medium: float=-0.002,
        drainage_threshold_high: float=0.003,
        drainage_threshold_medium: float=0.003,
    ):
    
    recharge = load_recharge(x, y, path_names)
    moving_mean_recharge = recharge['recharge (m/d)'].rolling(7, min_periods=1).mean().to_frame('moving mean recharge (m/d)')
    
    pssi_stage = pd.DataFrame(index = recharge.index)   
    if scenario == 'high':
        pssi_stage['scenario'] = np.where(moving_mean_recharge < infiltration_threshold_high, (surface_level-0.1),
                              np.where(moving_mean_recharge > drainage_threshold_high, (surface_level-0.5), summer_stage))
    elif scenario == 'medium':
        pssi_stage['scenario'] = np.where(moving_mean_recharge < infiltration_threshold_medium, (surface_level-0.1),
                              np.where(moving_mean_recharge > drainage_threshold_medium, (surface_level-0.5), summer_stage))
    else:
        ValueError('Not a valid "scenario": please define either "high" or "medium".')
    pssi_stage['scenario'] = pssi_stage['scenario'].where(~pssi_stage.index.month.isin([11,12,1,2]),(surface_level-0.5))    

    return pssi_stage[['scenario']]


def calculate_stats(surface_level, groundwater_level):
    relative_groundwater = surface_level - groundwater_level #surface level in mNAP and groundwaterleven in mNAP
    summer_mask = relative_groundwater['time.month'].isin([4,5,6,7,8,9])
    summer3_mask = relative_groundwater['time.month'].isin([6,7,8])
    winter_mask = relative_groundwater['time.month'].isin([10,11,12,1,2,3])
    
    summer_median = relative_groundwater.where(summer_mask, drop=True).mean(['time','x']).item()
    summer3_median = relative_groundwater.where(summer3_mask, drop=True).mean(['time','x']).item()
    winter_median = relative_groundwater.where(winter_mask, drop=True).mean(['time','x']).item()
    year_median = relative_groundwater.mean(['time','x']).item()
    
    rlg = relative_groundwater.quantile(0.9, 'time').mean(['x']).item()   
    rhg = relative_groundwater.quantile(0.1, 'time').mean(['x']).item()   
    
    return summer_median, summer3_median, winter_median, year_median, rlg, rhg #m - mv


#%%

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
            modflow_executable=to_linux(r'P:/11207812-somers-uitvoering/software/mf6_compiled_rm_v3'), #linux run
            #modflow_executable=to_linux(r'P:/11207812-somers-uitvoering/software/mf6.exe'), #local test
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
        output_10_perc = output['emission_be'].sum(['time', 'depth', 'x']).quantile(0.1, 'runs').item()
        output_90_perc = output['emission_be'].sum(['time', 'depth', 'x']).quantile(0.9, 'runs').item()
        output_min = output['emission_min'].sum(['time','depth','x']).median('runs').item()
        output_max = output['emission_max'].sum(['time','depth','x']).median('runs').item()
        
        min_summer_gwl,min_summer3_gwl,min_winter_gwl,min_year_gwl,min_rlg_gwl,min_rhg_gwl = calculate_stats(float(surface_level), output['phreatic_head'].min('runs'))
        mean_summer_gwl,mean_summer3_gwl,mean_winter_gwl,mean_year_gwl,mean_rlg_gwl,mean_rhg_gwl = calculate_stats(float(surface_level), output['phreatic_head'].median('runs'))
        max_summer_gwl,max_summer3_gwl,max_winter_gwl,max_year_gwl,max_rlg_gwl,max_rhg_gwl = calculate_stats(float(surface_level), output['phreatic_head'].max('runs'))

        nb_runs = len(output.runs)
        
        return (output_min, output_median, output_max, output_10_perc, output_90_perc,
                min_summer_gwl,min_summer3_gwl,min_winter_gwl,min_year_gwl,min_rlg_gwl,min_rhg_gwl,
                mean_summer_gwl,mean_summer3_gwl,mean_winter_gwl,mean_year_gwl,mean_rlg_gwl,mean_rhg_gwl,
                max_summer_gwl,max_summer3_gwl,max_winter_gwl,max_year_gwl,max_rlg_gwl,max_rhg_gwl, nb_runs)
        
    except (psycopg2.errors.UniqueViolation, ValueError, FloatingPointError, KeyError, FileNotFoundError, IndexError, TypeError, np.core._exceptions._ArrayMemoryError, PermissionError) as error:
        
        with open(to_linux(r"P:/11207812-somers-uitvoering/monitoring_2023/shp_files/error_lst_ref.txt"), "a") as file:
            file.write('parcel: '+str(parcel_id)+ "; " + str(error) + "\n")
            
        return (None, None, None, None, None, None,
                None, None, None, None, None, None, 
                None, None, None, None, None, None, 
                None, None, None, None, None, None)
              
#%%
def run_somers_ssi(parcel_id, surface_level, parcel_width, winter_stage, summer_stage, x, y, ssi_distance, ssi_depth, archetype): 
    parameter_combs = pd.read_csv(to_linux(r'P:/11207812-somers-ontwikkeling/somers_v2_development/data/3-input/pp2d_phreatic_head/parameters/definitive parameters/parameters_ssi_sv2.csv'), index_col = 0)
    try:
        output = somers.run_somers(
            ## model settings
            name = str(parcel_id),
            path_names = path_names,
            input_dir=to_linux(r'P:/11207812-somers-uitvoering/monitoring_2023/run_dir/3-input'),
            output_dir=to_linux(r'P:/11207812-somers-uitvoering/monitoring_2023/run_dir/4-output'),
            working_dir=to_linux(r'P:/11207812-somers-uitvoering/monitoring_2023/run_dir/2-interim'),
            modflow_executable=to_linux(r'P:/11207812-somers-uitvoering/software/mf6_compiled_rm_v3'), #linux
            #modflow_executable=to_linux(r'P:/11207812-somers-uitvoering/software/mf6.exe'), #local test
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
            parcel_type = 'ssi',
            x=x,
            y=y,
            surface_level = float(surface_level),
            parcel_width = float(parcel_width),
            start_date=f'01-01-{find_weerjaar(x, y)}',
            end_date=f'31-12-{find_weerjaar(x, y)}',
        
            ## ditch parameters
            ditch_depth = find_ditch_depth(surface_level, np.min([winter_stage, summer_stage])),
            winter_stage = float(winter_stage),
            summer_stage = float(summer_stage),
                    
            ## ssi-system parameters
            drain_distance = float(ssi_distance),
            drain_depth = float(ssi_depth), 
    
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
        output_10_perc = output['emission_be'].sum(['time', 'depth', 'x']).quantile(0.1, 'runs').item()
        output_90_perc = output['emission_be'].sum(['time', 'depth', 'x']).quantile(0.9, 'runs').item()
        output_min = output['emission_min'].sum(['time','depth','x']).median('runs').item()
        output_max = output['emission_max'].sum(['time','depth','x']).median('runs').item()
        
        min_summer_gwl,min_summer3_gwl,min_winter_gwl,min_year_gwl,min_rlg_gwl,min_rhg_gwl = calculate_stats(float(surface_level), output['phreatic_head'].min('runs'))
        mean_summer_gwl,mean_summer3_gwl,mean_winter_gwl,mean_year_gwl,mean_rlg_gwl,mean_rhg_gwl = calculate_stats(float(surface_level), output['phreatic_head'].median('runs'))
        max_summer_gwl,max_summer3_gwl,max_winter_gwl,max_year_gwl,max_rlg_gwl,max_rhg_gwl = calculate_stats(float(surface_level), output['phreatic_head'].max('runs'))

        nb_runs = len(output.runs)
        
        return (output_min, output_median, output_max, output_10_perc, output_90_perc,
                min_summer_gwl,min_summer3_gwl,min_winter_gwl,min_year_gwl,min_rlg_gwl,min_rhg_gwl,
                mean_summer_gwl,mean_summer3_gwl,mean_winter_gwl,mean_year_gwl,mean_rlg_gwl,mean_rhg_gwl,
                max_summer_gwl,max_summer3_gwl,max_winter_gwl,max_year_gwl,max_rlg_gwl,max_rhg_gwl, nb_runs)
        
    except (psycopg2.errors.UniqueViolation, ValueError, FloatingPointError, KeyError, FileNotFoundError, IndexError, TypeError, np.core._exceptions._ArrayMemoryError, PermissionError) as error:
        
        with open(to_linux(r"P:/11207812-somers-uitvoering/monitoring_2023/shp_files/error_lst_ssi.txt"), "a") as file:
            file.write('parcel: '+str(parcel_id)+ "; " + str(error) + "\n")
            
        return (None, None, None, None, None, None,
                None, None, None, None, None, None, 
                None, None, None, None, None, None, 
                None, None, None, None, None, None)
    
#%%
def run_somers_pssi(parcel_id, surface_level, summer_stage, winter_stage, parcel_width,x, y, pssi_distance, pssi_depth, archetype): 
    parameter_combs = pd.read_csv(to_linux(r'P:/11207812-somers-ontwikkeling/somers_v2_development/data/3-input/pp2d_phreatic_head/parameters/definitive parameters/parameters_ssi_sv2.csv'), index_col = 0)
    preset_pssi_stage = set_pssi_scenario(scenario ='medium',x=x,y=y,surface_level=float(surface_level),summer_stage=float(summer_stage),path_names=path_names)
    try:
        output = somers.run_somers(
            ## model settings
            name = str(parcel_id),
            path_names = path_names,
            input_dir=to_linux(r'P:/11207812-somers-uitvoering/monitoring_2023/run_dir/3-input'),
            output_dir=to_linux(r'P:/11207812-somers-uitvoering/monitoring_2023/run_dir/4-output'),
            working_dir=to_linux(r'P:/11207812-somers-uitvoering/monitoring_2023/run_dir/2-interim'),
            modflow_executable=to_linux(r'P:/11207812-somers-uitvoering/software/mf6_compiled_rm_v3'), #linux
            #modflow_executable=to_linux(r'P:/11207812-somers-uitvoering/software/mf6.exe'), #local test
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
            parcel_type = 'pssi',  
            x=x,
            y=y,
            surface_level = float(surface_level),
            parcel_width = float(parcel_width),
            start_date=f'01-01-{find_weerjaar(x, y)}',
            end_date=f'31-12-{find_weerjaar(x, y)}',
        
            ## ditch parameters
            ditch_depth=find_ditch_depth(surface_level, np.min([winter_stage, summer_stage])),
            winter_stage = float(winter_stage),
            summer_stage = float(summer_stage),
                    
            ## ssi-system parameters
            drain_distance = float(pssi_distance),
            drain_depth = float(pssi_depth), 
            pssi_winter_stage = float(winter_stage), #scenario aanpak toevoegen
            pssi_summer_stage = float(summer_stage),
            
            ## preset values
            preset_soilcode = archetype, 
            preset_pssi_stage = preset_pssi_stage,
            
            ## output 
            save_input = False,
            save_modflow_runs = False,
            save_phreatic_head = False,
            save_soil_moisture = False,
            save_soil_temperature = False,
            save_emission= False,
            
            )
        
        output_median = output['emission_be'].sum(['time','depth','x']).median('runs').item()
        output_10_perc = output['emission_be'].sum(['time', 'depth', 'x']).quantile(0.1, 'runs').item()
        output_90_perc = output['emission_be'].sum(['time', 'depth', 'x']).quantile(0.9, 'runs').item()
        output_min = output['emission_min'].sum(['time','depth','x']).median('runs').item()
        output_max = output['emission_max'].sum(['time','depth','x']).median('runs').item()
        
        min_summer_gwl,min_summer3_gwl,min_winter_gwl,min_year_gwl,min_rlg_gwl,min_rhg_gwl = calculate_stats(float(surface_level), output['phreatic_head'].min('runs'))
        mean_summer_gwl,mean_summer3_gwl,mean_winter_gwl,mean_year_gwl,mean_rlg_gwl,mean_rhg_gwl = calculate_stats(float(surface_level), output['phreatic_head'].median('runs'))
        max_summer_gwl,max_summer3_gwl,max_winter_gwl,max_year_gwl,max_rlg_gwl,max_rhg_gwl = calculate_stats(float(surface_level), output['phreatic_head'].max('runs'))

        nb_runs = len(output.runs)
        
        return (output_min, output_median, output_max, output_10_perc, output_90_perc,
                min_summer_gwl,min_summer3_gwl,min_winter_gwl,min_year_gwl,min_rlg_gwl,min_rhg_gwl,
                mean_summer_gwl,mean_summer3_gwl,mean_winter_gwl,mean_year_gwl,mean_rlg_gwl,mean_rhg_gwl,
                max_summer_gwl,max_summer3_gwl,max_winter_gwl,max_year_gwl,max_rlg_gwl,max_rhg_gwl, nb_runs)
        
    except (psycopg2.errors.UniqueViolation, ValueError, FloatingPointError, KeyError, FileNotFoundError, IndexError, TypeError, np.core._exceptions._ArrayMemoryError, PermissionError) as error:
        
        with open(to_linux(r"P:/11207812-somers-uitvoering/monitoring_2023/shp_files/error_lst_pssi.txt"), "a") as file:
            file.write('parcel: '+str(parcel_id)+ "; " + str(error) + "\n")
            
        return (None, None, None, None, None, None,
                None, None, None, None, None, None, 
                None, None, None, None, None, None, 
                None, None, None, None, None, None)
