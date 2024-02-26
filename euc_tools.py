import xarray as xr; import numpy as np; 
import pandas as pd; 
import thermo

def vertical_chunk( ds, zmax ):
    # Cut a vertical chunk of POP2 data across all z_t, z_w_top, and z_w_bot
    ds = ds.sel( z_t = slice( 0, zmax ) ); 
    Nz = len( ds['z_t'] ); # how many left
    ds = ds.isel( z_w_top = [jj for jj in range( Nz ) ] )
    ds = ds.isel( z_w_bot = [jj for jj in range( Nz ) ] )
    return ds 

# Tools to describe the EUC
def get_euc_area( ds ):
    ds = ds.sel( lat = slice( -3, 3 ) ); # equatorial
    ds = ds.sel( lon = slice( 125, 280 ) ); # pacific
    ds = vertical_chunk( ds, 500 ); # cut from surface to 500 m
    return ds 

def find_euc( ds ):
    # Take ocean dataset from pop2 and apply conditions to 
    # return a mask that is True in locations that are the EUC
    # Begin with spatial subsetting
    ds = get_euc_area( ds ); 
    # Now physical / ocean-dependent conditions
    is_eastward = ds['UVEL'] > 10; # in cm / s 
    below_ML = ds['z_t'] > ( ds['HMXL']/100 );
    # Combine masks
    euc_mask = is_eastward * below_ML
    return ds.where( euc_mask ) 
 
def get_euc_transport( ds, sum_over = ['lat','z_t'] ):
    ds = find_euc( ds ); # subselect EUC
    # Get numbers necessary for integrating transport
    dx, dy = thermo.get_dx_dy( ds ); 
    dz = thermo.get_dz( ds ) 
    # Now integrate
    transport = ( ds['UVEL'] / 100 * dy * dz / 1e6 ).sum( sum_over ).persist()
    return transport 

def find_euc_bottom( ds ):
    # FOR SIMPLICITY I THINK THIS WILL HAVE TO BECOME FIND EUC_CORE
    # Start by calculating EUC transport without summing over z_t
    transport = get_euc_transport( ds , sum_over = ['lat'] );
    # Already multiplied by dy and dz. 
    vert_int = ( transport * dz ).sum(['z_t']); # total transport across all depth
    # Now find depths above which 95% of transport is
    cum_trans = transport.cumsum( ['z_t'] ); 

def find_euc_core( ds ):
    # Get mean-lat of euc uvel, and find depth of maximum velocity
    ds = find_euc( ds ); 
    current = ds['UVEL'].mean('lat'); 
    ind_z = current.argmax( 'z_t' ); # depth index of greatest velocity
    return ind_z

def euc_plane_view( ds ):
    # Start by finding EUC core and its density
    core_z = find_euc_core( ds ); 
    # More lenient spatial subsetting
    ds = ds.sel( lon = slice( 135, 290 ) ).sel( lat = slice( -25, 25 ) ); 
    ds = vertical_chunk( ds , 500 ); 
    # Now find density of EUC core
    core_dens = ds['RHO'].sel( lat = slice(-3,3) ).mean('lat').isel( z_t = core_z ); # density of EUC core
    # find regions that can feed the EUC core
    core_feeder = ds.where( np.abs( ds['RHO'] - core_dens ) < 0.001 ).where( ds['z_t'] > ds['HMXL'] * 100 ); 
    dz = thermo.get_dz( ds );
    # Get plane view
    core_feeder = core_feeder.weighted( dz ).mean('z_t'); # needs further testing
    return core_feeder
        
def KE_budget( ds ):
    # Evaluate terms in the kinetic energy budget of the EUC
    #ds = get_euc_area( ds )
    # Basic ingredients
    u = ( ds['UVEL']/100 ).persist()
    dudz = - u.differentiate('z_t').persist() ;
    rho = 1024; # kg m^-3
    # Calculate viscous and vertical-advective components of KE change
    vert = - ( u * ds['WVEL'].values/100 * dudz * rho ).persist()
    visc = ( - rho * u * ( ds['VVC'].values * 1e-4  * dudz ).differentiate('z_t') ).persist()
    return  vert, visc 
    


