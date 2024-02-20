import xarray as xr; import numpy as np; 
import pandas as pd; 
import thermo

# Tools to describe the EUC
def get_euc_area( ds ):
    ds = ds.sel( lat = slice( -3, 3 ) ); # equatorial
    ds = ds.sel( lon = slice( 125, 280 ) ); # pacific
    ds = ds.sel( z_t = slice( 0, 500 ) );
    # To continue vertical subsetting use:
    Nz = len( ds['z_t'] ); 
    ds = ds.isel( z_w_top = [ jj for jj in range( Nz ) ] )
    ds = ds.isel( z_w_bot = [ jj for jj in range( Nz ) ] )
    return ds 

def find_euc( ds ):
    # Take ocean dataset from pop2 and apply conditions to 
    # return a mask that is True in locations that are the EUC
    
    # Begin with spatial subsetting
    ds = get_euc_area( ds ); 
    # Now physical / ocean-dependent conditions
    is_eastward = ds['UVEL'] > 5; # in cm / s 
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
    # Start by calculating EUC transport without summing over z_t
    transport = get_euc_transport( ds , sum_over = ['lat'] ); 
    vert_int = transport.sum(['z_t']);
    # Now use cumsum to find depths above which 95% of transport is
    cum_trans = transport.cumsum( ['z_t'] ); 
    pass

