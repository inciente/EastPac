import xarray as xr; import pandas as pd; 
import gsw; import numpy as np; 

'''
Functions to deal with the thermodynamics of ocean and atmospheric components
of model simulations. 

Might have to refactor and replace coordinate references for flexible ones using model instances from intercomparison.

'''

def find_z_temp( temp_xr, T = 20 ):
    # Find depth at which temperature is closest to T
    temp_xr = temp_xr.fillna( -300 ); # replace nans
    depth_index = np.abs( temp_xr - T ).argmin( dim = 'z_t' )
    # Extract corresponding depth values
    z_temp = ( temp_xr['z_t'] ).isel( z_t = depth_index )
    # Mask out locations where value is at surface (helps with nans)
    z_temp = xr.where( z_temp > 6 , z_temp, np.nan )
    return z_temp

def area_weights( xr_obj ):
    x = 'lon'; y = 'lat';
    # Create weights to compensate with latitude
    weights = np.cos( np.deg2rad( xr_obj[ y ] ) );
    weights.name = 'weights'
    return weights

def get_dx_dy( xr_obj ):
    dx = np.abs( np.gradient( xr_obj['lon'] ) );
    dx = xr.DataArray( data = dx, coords = {'lon':xr_obj['lon'] } )
    dy = np.abs( np.gradient( xr_obj['lat'] ) ); 
    dy = xr.DataArray( data = dy, coords = {'lat':xr_obj['lat'] } )
    return dx, dy

def get_dz( xr_obj , source = 'ocn' ):
    dz = np.abs( xr_obj['z_w_top'].values - xr_obj['z_w_bot'].values )
    dz = xr.DataArray( data = dz, coords = {'z_t':xr_obj['z_t']} )
    return dz 

def heat_per_cell( ocn_xr ):
    # Amount of heat in each grid cell, relative to T = 0 Celsius
    lay_th = get_dz( ocn_xr )
    c_p = 3.9e3; 
    dx, dy = get_dx_dy( ocn_xr ); # get grid info
    # calculate heat density
    heat_in_grid = c_p * ocn_xr['RHO'] * ocn_xr['TEMP'] 
    heat_in_grid = heat_in_grid * lay_th * dx * dy ;# times meters * degrees ** 2 
    heat_in_grid = heat_in_grid * 110e3 ** 2 * area_weights( ocn_xr ); # to meters
    # Still need to integrate 
    return heat_in_grid 







