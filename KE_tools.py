import sys; 
import matplotlib.pyplot as plt; 
import numpy as np; 
import pandas as pd; import xarray as xr; 
sys.path.append('/home/noelgb/repositories/EastPac/')
import intercomparison, thermo
sys.path.append('/home/noelgb/notebooks/config_files')
from open_jsons import *
# --------------

''' 
Functions to compute the energy budget of the upper Pacific Ocean
'''

def meridional_KE_flux( pop_ds, integrate = True ):
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

def meridional_pwork( pop_ds, ref_rho, ref_ssh, integrate = True  ):
    # Compute pressure the meridional pressure work
    vel = pop_ds[ 'VVEL' ] / 100
    # specify reference density and ssh 
    p_anom = thermo.pressure( pop_ds , ref_rho = ref_rho , 
               ssh_val = ref_ssh )
    pwork = vel * p_anom
    if integrate:
        pwork = xz_integral( pwork, pop_ds )
    return pwork

def meridional_APE_flux( pop_ds, ref_rho, integrate = True ):
    # compute APE and its meridional flux
    rho_prime = pop_ds['RHO'] * 1000 - ref_rho; 
    N2_ref = 9.81 / 1025 * ref_rho.differentiate( 'z_t' )
    APE = 9.81 ** 2 * rho_prime ** 2 
    APE = APE / ( 2 * 1025 * N2_ref )
    APE_flux = pop_ds['VVEL'] / 100 * APE
    
    if integrate:
        APE_flux = xz_integral( APE_flux, pop_ds )
    return APE_flux

def mixing_work( pop_ds, integrate = True ):
    # volume integral of tpower. assuming there's been some
    # spatial subsetting (or application of where)
    power = pop_ds['TPOWER'] * 1e-7 * 1e6 # to J / m3 / s
    if integrate: 
        power = xyz_integral( power.rename( {'z_w_bot':'z_t'} ),
                               pop_ds )
    return power

def wind_work( pop_ds, integrate = True ):
    # wind work on the ocean surface
    power = pop_ds['TAUX'] * pop_ds['UVEL'].isel( z_t = 0 )
    power += pop_ds['TAUY'] * pop_ds['VVEL'].isel( z_t = 0 )
    # now convert to W m-2
    power = power * 1e-7 * 1e4
    if integrate:
        power = xy_integral( power, pop_ds )
    return power  
 
def upwelling_work( pop_ds, ref_rho , integrate = True ):
    # energy cost of upwelling
    rho = pop_ds['RHO'] * 1000 - ref_rho
    power = rho * 9.81 * pop_ds['WVEL'].values / 100
    if integrate:
        power = xyz_integral( power , pop_ds )
    return power

# -------- functions for spatial subsetting and section management

def basin_mask():
    # Load the WOA basin mask at 100 m depth
    path = '/home/noelgb/data/ocean_mask/WOA_basins_100m.nc'
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
    mask = mask['basins'].interp( lat = xr_obj['lat'] ).interp( \
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
    ds = pacific_only( ds );
    ds = ds.isel( time = [0,12,24] ); # for speed up in testing
    # get reference density and ssh
    ds['rho_ref'] = thermo.reference_density( ds )
    ds['rho_ref'] = ds['rho_ref'].interp( z_t = ds['z_t'] , method = 'linear', 
                                         kwargs = {'fill_value': 'extrapolate' } )
    ds['mean_ssh'] = ds['SSH'].mean(['lon','lat','time']).persist()
    return ds 

def all_merid_fluxes( pop_ds ):
    # compute all meridional fluxes from functions in this module
    fluxes = xr.Dataset()
    fluxes['adv'] = meridional_KE_flux( pop_ds , integrate = True );
    fluxes['pwork'] = meridional_pwork( pop_ds , 
                         ref_rho = pop_ds['rho_ref'], 
                         ref_ssh = pop_ds['mean_ssh'],
                         integrate = True )
    fluxes['APE_adv'] = meridional_APE_flux( pop_ds , 
                         ref_rho = pop_ds['rho_ref'], 
                         integrate = True )
    return fluxes 

def all_conversions( pop_ds ):
    # compute all energy conversions involving area/volume integrals
    # each function returns a time series
    conversions = xr.Dataset()
    conversions['upwelling'] = upwelling_work( pop_ds , 
                                    ref_rho = pop_ds['rho_ref'],
                                    integrate = True )
    conversions['windwork'] = wind_work( pop_ds , integrate = True )
    conversions['mixing'] = mixing_work( pop_ds, integrate = True )
    return conversions

def eval_KE( model ):
    ds = prepare_model_for_energy( model )
    latitudes = [-35, -15, -3, 3, 15, 35 ]
    # Get all the components of the energy budget
    
    # ---------- begin with meridional fluxes
    ds_merid = ds.sel( lat = latitudes , method = 'nearest' );
    fluxes = all_merid_fluxes( ds_merid ); 
    
    # ---------- now slice sections where area or volume ints involved
    area_quants = []
    for jj in range( len( latitudes ) - 1 ):
        lat_range = slice( latitudes[jj] , latitudes[jj+1] )
        ds_here = ds.sel( lat = lat_range ); # this will be integrated
        conversions = all_conversions( ds_here )
        conversions = conversions.assign_coords( { 'lat': [latitudes[jj]] } )
        area_quants.append( conversions )
    # concatenate all area and volume integrals
    conversions = xr.concat( area_quants, dim = 'lat' )
    return fluxes, conversions 


    

    


