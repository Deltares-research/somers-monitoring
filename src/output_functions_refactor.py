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

#%% SOMERS functions

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
    weerregios = gpd.read_file(r'P:/11207812-somers-ontwikkeling/somers_v2.1_input/data/meteo/weather_regions.shp')
    weerregio = weerregios[weerregios.contains(point)]['weather_rg'].item()
    if weerregio == 'southwest':
        return '2015'
    elif weerregio == 'northeast':
        return '2013'

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