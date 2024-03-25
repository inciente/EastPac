import sys; 
import matplotlib.pyplot as plt; 
import numpy as np; 
import pandas as pd; import xarray as xr; 
sys.path.append('/home/noelgb/repositories/EastPac/')
import intercomparison, thermo
sys.path.append('/home/noelgb/config_files')
from open_jsons import *
# --------------

''' 
Functions to compute the energy budget of the upper Pacific Ocean
'''

def KE_flux( pop_ds, integrate = True ):
    # computes KE advective flux in given direction
    # pop_ds is xr.dataset with pop2 output
    # direction determines advective variable
    direction = 'VVEL'
    KE = 1025 / 2 * ( pop_ds['UVEL'] ** 2 + pop_ds['VVEL'] ** 2 \
           + pop_ds['WVEL'].values ** 2 ) * pop_ds[ direction ] / 1e6
    # division over 1e6 adjusts to m/s
    if integrate:
        KE = xz_integral( KE, pop_ds )
    return KE

def merid_pwork( pop_ds, ref_rho, ref_ssh, integrate = True  ):
    # Compute pressure the meridional pressure work
    vel = pop_ds[ 'VVEL' ] / 100
    # specify reference density and ssh 
    p_anom = thermo.pressure( pop_ds , ref_rho = ref_rho , 
               ssh_val = ref_ssh )
    pwork = vel * p_anom
    if integrate:
        pwork = xz_integral( pwork, pop_ds )
    return pwork

def mixing_work( pop_ds ):
    # volume integral of tpower. assuming there's been some
    # spatial subsetting (or application of where)
    power = pop_ds['TPOWER'] * 1e-7 * 1e6 # to J / m3 / s
    if integrate: 
        power = xyz_integral( power.rename( {'z_w_bot':'z_t'} ),
                               pop_ds )
    return power

def wwork( pop_ds, integrate = True ):
    # wind work on the ocean surface
    power = pop_ds['TAUX'] * pop_ds['UVEL'].isel( z_t = 0 )
    power += pop_ds['TAUY'] * pop_ds['VVEL'].isel( z_t = 0 )
    # now convert to W m-2
    power = power * 1e-7 * 1e4
    if integrate:
        power = xy_integral( power, pop_ds )
    return power  
 
def upwork( pop_ds, ref_rho , integrate = True ):
    # energy cost of upwelling
    rho = pop_ds['RHO'] * 1000 - ref_rho
    power = rho * 9.81 * pop_ds['WVEL'].rename( {'z_w_top':'z_t'} )
    if integrate:
        power = xyz_integral( power )
    return power

# -------- functions for spatial subsetting and section management

def basin_mask():
    # Load the WOA basin mask at 100 m depth
    path = '/home/noelgb/data/ocean_masks/WOA_basins_100m.nc'
    mask = xr.open_dataset( path ); 
    # change longitude from -180 to 0-360
    nulon = mask['lon'].values;
    negative = nulon < 0; 
    nulon[negative] = nulon[negative] + 360;
    mask['lon'] = nulon;
    # sort 
    mask = mask.sortby('lon')
    return mask

def pacific_only( xr_obj ):
    # Apply WOA mask to object, leaving only data inside the Pacific
    mask = basin_mask(); 
    mask = mask.interp( lat = xr_obj['lat'] ).interp( \
                        lon = xr_obj['lon'] )
    return xr_obj.where( mask == 2 )

def xy_integral( xr_obj, parent_ds ):
    # return area integral (lon, lat)
    data = ( xr_obj * parent_ds['dx'] * parent_ds['dy'] )
    data = data.sum(['lon','lat'])
    return data

def xz_integral( xr_obj, parent_ds ):
    # integrate over lon, depth
    data = ( xr_obj * parent_ds['dx'] * parent_ds['dz'] )
    data = data.sum( ['lon','z_t'] )
    return data

def z_integral( xr_obj, parent_ds ):
    data = ( xr_obj * parent_ds['dz'] ).sum('z_t')
    return data

def xyz_integral( xr_obj , parent_ds ):
    data = xy_integral( xr_obj, parent_ds )
    data = z_integral( data, parent_ds )
    return data 

'''
LOOKS LIKE ITS GOING TO BE EASIER TO DO WITHOUT THIS CLASS
'''
#class lat_line:
#    
#    def __init__( self , lat, xlims, ds ):
#        self.lat = lat; 
#        self.xlims = xlims;
#        #self.ds = self.eval_xr( ds ); # parent xr to inherit info
#
#    def eval_xr( self, xr_obj ):
#        # Apply spatial subsetting 
#        xr_obj = xr_obj.sel( lon = self.xlims ).sel( lat = self.lat,
#                  method = 'nearest' )
#        return xr_obj 
#
#    def area_integral( self, xr_obj, parent_ds ):
#        xr_obj = self.eval_xr( xr_obj )
#        # integrate over depth and longitude
#        data = ( xr_obj * self.ds['dx'] * self.ds['dy'] ).sum( \
#                 ['lon','z_t'] )


''' 
NOW ACTUALLY LAY OUT THE COMPUTATION WE'RE INTERESTED IN
'''

def prepare_model_for_energy( model ):
    ds = pull_combined_json( model )
    latitudes = [-35,-15,-3,3,15,35];
    # get metrics and spatial subsetting
    ds = ds.sel( lat = slice( -40, 40 ) );
    ds = ds.sel( lon = slice( 100, 295 ) );
    ds = pacific_only( ds )
    # get reference density and ssh
    ds['rho_ref'] = thermo.reference_density( ds )
    ds['mean_ssh'] = ds['SSH'].mean(['lon','lat','time']).persist()
    return ds 

def eval_KE( model ):
    ds = prepare_model_for_energy( model )
    latitudes = [-30, -15, -3, 3, 15, 30 ]
    # Get all the components of the energy budget
    pass

    


