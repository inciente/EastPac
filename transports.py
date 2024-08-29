# Standard Python packages
import xarray as xr; import numpy as np; 
import pandas as pd; import sys, os, time;
import shapely as shp

# Other packages within EastPac repo
import thermo; 

''' 
Suite of tools to calculate energy transports in POP2 and CAM6 model output.

'''

# --------------- OCEANIC ROUTINES -------------------   

def ocean_diabatic_heating( ds ): 
    # Input is xr.Dataset from POP2
    dz = thermo.get_dz( ds );
    
    # Add up temperature tendencies and advective contribution
    diabatic = 1024 * 3.9e3 * ( ds['TEND_TEMP'] - ds['ADV_3D_TEMP'] )
    diabatic = ( ds['dz'] * diabatic ).sum( 'z_t' )
    return diabatic

def ocn_zonal_integral( ds, xr_obj ):
    # Return zonal integral of ocean variables.
    # use ULONG for greater precision at high latitudes
    dx = 110e3 * ds['ULONG'].isel( time = 0 ).differentiate( 'lon' ) * thermo.area_weights( ds )
    integral = ( dx * xr_obj ).sum( 'lon' )
    return integral

# ------------ ATMOSPHERIC ROUTINES ------------------

def atm_vertical_integral( ds, xr_obj ):
    # Integrate xr_obj with respect to pressure, dividing by g.
    # grid information is stored in ds
    layer_thickness = xr.DataArray( data = np.diff( ds['ilev'].values ),                         dims = ('lev') , coords = {'lev':ds['lev'] } )
    vert_integ = ( layer_thickness * xr_obj ).sum( 'lev' ) / 9.81
    return vert_integ 

def ocn_zonal_integral( ds, xr_obj ):
    dx, dy = thermo.get_dx_dy( ds )
    

def atm_OLR( ds , integrate = True ):
    # Input is xr.Dataset from CAM6
    # get upwelling longwave radiation at top of model
    olr = ds['FLUT']
    if integrate:
        dx, dy = thermo.get_dx_dy( ds );
        olr = ( dx * olr ).sum('lon')
    return olr 

def atm_net_solar( ds, integrate = True ):
    # Input is xr.Dataset from CAM6
    # get net solar flux at top of atmosphere
    toa = ds['FSNTOA']
    if integrate:
        dx, dy = thermo.get_dx_dy( ds );
        toa = ( dx * toa ).sum( 'lon' )
    return toa
    
def atm_mse_transport( ds, integrate = True ):
    # Input is xr.Dataset from CAM6
    # return meridional transport of moist static energy

    # necessary constants
    Lv = metpy.constants.water_heat_vaporization.magnitude
    Cpd = metpy.constants.dry_air_spec_heat_press.magnitude
    g = 9.81

    # estimate mse transport
    htrans = g * ds['Z3'] + Cpd * ds['VT'] + Lv * ds['VQ']
    if integrate: 
        htrans = atm_vertical_integral( ds, htrans )
    return htrans 
    







