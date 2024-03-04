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
    # Get xy grid with dx and dy for xarray object in meteres
    # use for diffrentiation or for computing area integrals
    xx, yy = np.meshgrid( xr_obj['lon'], xr_obj['lat'] );
    dims = ('lat','lon'); 
    coords = {'lon':xr_obj['lon'], 'lat':xr_obj['lat']}
    # Translate to xarray to maintain dimensions
    xx = xr.DataArray( xx, dims = dims, coords = coords )
    yy = xr.DataArray( yy, dims = dims, coords = coords );
    try:
        dx = 110e3 * xx.differentiate('lon') * area_weights(xr_obj)
    except: 
        dx = None
    try:    
        dy = 110e3 * yy.differentiate('lat')
    except:
        dy = None; # in case object has no lat 
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
    #heat_in_grid = heat_in_grid * 110e3 ** 2 * area_weights( ocn_xr ); # to meters
    # Still need to integrate 
    return heat_in_grid 

def pacific_gradient( xr_obj ):
    # Use xr_obj to compute gradient between west and east pacific
    EP_lon = slice(210, 270); WP_lon = slice( 140, 155 )
    xr_obj = xr_obj.sel( lat = slice( -5, 5 ) ).mean( dim = 'lat' );
    # Get averages for each sector
    WP_xr = xr_obj.sel( lon = WP_lon ).mean( dim = 'lon' )
    EP_xr = xr_obj.sel( lon = EP_lon ).mean( dim = 'lon' );
    return (EP_xr - WP_xr)

def pressure( ds , add_atm = None ):
    g = 9.81; rho = ds['RHO'] * 1000; # fix units
    # Take output from POP2 (ds) and compute pressure everywhere
    psurf = ds['SSH'] / 100 * g * rho.isel( z_t = 0 ); 
    # Option to add atmospheric pressure
    if add_atm is not None:
        psurf = psurf + add_atm; # this assumes interpolation of PSL to ocean grid
    # Get dz to integrate density
    dz = get_dz( ds ); 
    # Compute weight of each model cell
    weight = ( rho * g * dz )
    # Now write the pressure sum
    pressure = psurf + weight.cumsum( 'z_t' ).shift( z_t = 1 ,
                fill_value = 0 ) + weight / 2; 
    # cumsum.shift integrates weight of cells above only
    # weight/2 accoutns for weight of current cell 
    pressure.attrs['units'] = 'Pa [N m-2]'
    pressure.attrs['description'] = 'Ocean pressure at depth z_t'
    return pressure 





