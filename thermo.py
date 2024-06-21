import xarray as xr; import pandas as pd; 
import gsw; import numpy as np; 
from shapely import Polygon, Point


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
    z_temp = ( temp_xr['z_t'] ).isel( z_t = depth_index.compute() )
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

def pressure( ds , add_atm = None, ref_rho = None, 
           ssh_val = None):
    g = 9.81; rho = ds['RHO'] * 1000; # fix units
    
    if ref_rho is not None:
        # in case we want to compute pressure anomaly
        rho = rho - ref_rho; 
    # Take output from POP2 (ds) and compute pressure everywhere 
    if ssh_val is not None:
        ds['SSH'] = ds['SSH'] - ssh_val;
    psurf = ds['SSH'] / 100 * g * 1022;
    
    #if ref_rho is not None:
    #    if ssh_val is not None:
    #        psurf = psurf - ssh_val 
    #    # a;sp get anomaly of ssh
    #    psurf = psurf - psurf.mean( ['lon','lat'] ) 
    # Option to add atmospheric pressure

    if add_atm is not None:
        psurf = psurf + add_atm; # this assumes interpolation of PSL to ocean grid
    # Get dz to integrate density
    try:
        dz = ds['dz']
    except:
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

def polygon_mask( xr_obj , poly, dims = ['lon','lat'] ):
    # Check what points in the xr_obj grid are within the polygon
    # and return map of booleans.
    mx, my = np.meshgrid( xr_obj[ dims[0] ] , xr_obj[ dims[1] ] ); 
    ar_shape = mx.shape; # to reshape final product
    mx = mx.flatten(); my = my.flatten();

    # create points from mx, my coords and compare to polygon
    points = [ Point( ( mx[jj], my[jj] ) ) for jj in range( len( mx ) ) ];
    inside = np.array( [ pp.within( poly ) for pp in points ] ); # booleans
    # --------- give inside the right dimensions and store as xarray
    inside = xr.DataArray( data = inside.reshape( ar_shape ),
                           coords = { dims[1] : xr_obj[ dims[1] ], dims[0] : xr_obj[ dims[0] ] } )
    return inside

def add_measures_of_volume( ds ):
    ds['dz'] = get_dz( ds ); 
    ds['dx'], ds['dy'] = get_dx_dy( ds )
    return ds 

def mass_histograms( ds ):
    # Return histograms of density and height, weighted by cell volume
    stacker = {'space' : [ dim for dim in ds.dims if dim in ['z_t','lat','lon'] ] };

    # reshape density, and heights available
    #rho_to_order = 1000 * ds['RHO'].stack( stacker );
    sigma_to_order = 1000 + gsw.sigma0( ds['SALT'], ds['TEMP'] )
    sigma_to_order = sigma_to_order.stack( stacker );
    z_available = (ds['z_t'] * ds['RHO']/ds['RHO']).stack( stacker)
    mask = ~ np.isnan( sigma_to_order ).compute(); # remove nans

    # compute cell volume to weight space and mass available
    cell_volume = (ds['dx'] * ds['dy'] * ds['dz']).stack( stacker )
    
    # define bins for histograms
    rho_bins = np.arange( 1012.5, 1065, 0.25 ); 
    z_bins = xr.concat( [ ds['z_t'] - 1 , ds['z_t'].isel( \
                   z_t = -1 ) + 1 ], dim = 'z_t' )

    # compute histograms
    #rho_hist = np.histogram( rho_to_order[mask] , bins = rho_bins,
    #               weights = cell_volume[mask] )
    
    sigma_hist = np.histogram( sigma_to_order[mask] , bins = rho_bins , 
                  weights = cell_volume[mask] );

    z_hist = np.histogram( z_available[mask] , bins = z_bins,
                   weights = cell_volume[mask] );
    
    #return rho_hist, z_hist
    return sigma_hist, z_hist

def rho_from_sigma( sigma, pressure ):
    # derived from least square fit. meant to be used in reference_density
    # sigma must be 1020-ish (not 20-ish)
    # pressure must be in dbar (or meters)
    rho = sigma * 1.052076 + pressure * 0.004589025 - 53.3224
    return rho 

def divide_surface_layer( ds, z_hist ):
    '''
     take z output from mass_histograms, and add a super-thin
     layer in the surface. this helps allocate super-light water
    '''
    zbins = z_hist[1]; # bins
    z_vol = z_hist[0]; # volume within each bin
    surf_vol = z_vol[0]; 
    # volume in surface layer, which is defined by z limits [ zbins[0], zbins[1] ] 
    surf_thickness = ds['dz'][0].values; # thickness of surface layer
    nu_bins = np.array( [ 0, 0.0001, 0.001, 0.005, 0.01 , 0.02 , 
            0.05, 0.1, 0.25, 0.5, 0.8, 1 , zbins[0] ] ); # new zlims to add before zbins
    nu_dz = np.diff( nu_bins ); # thickness of new layers (this well completely replace surface layer)
    # ---- 
    volume_left = ( surf_thickness - nu_bins[-1] ) / surf_thickness * surf_vol; # needs to be further divided in two 
    z_vol[0] = volume_left; # update volume of surface layer that shrunk
    # -------- now allocate lost volume to new layers
    new_volume = surf_vol - volume_left; # to be allocated across new layers
    allocated_volumes = nu_dz / ( nu_bins[-1]-nu_bins[0] )  * new_volume
    # update z_hist
    z_hist[1] = np.concatenate( [ nu_bins[:-1] , zbins ] ); # bin edges
    z_hist[0] = np.concatenate( [ allocated_volumes, z_vol ] );
    
    return z_hist


def reference_density( ds , subset = None ): 
    # Rearrange density to find profile with minimum PE
    ds = add_measures_of_volume( ds )    
    if subset is not None:
        ds = subset( ds ); # will get warning if not callable
    ds = ds[ ['RHO','TEMP','SALT','dz','dy','dx'] ];
    try:
        ds = ds.mean('time'); # can't do for each step
    except: 
        pass 
    # compute pdfs of rho and z weighted by volume
    # current version does pdf of sigma... in development 
    #rhoh, zh = mass_histograms( ds )
    sigmah, zh = mass_histograms( ds );
    zh = divide_surface_layer( ds, list( zh ) ); # update
    # cross-correlate cumulative distributions of height and rho
    # --- kill tail to avoid super light water taking up surface
    cumulative = np.linspace( 0, sum(zh[0]) ,
                     1500 ); 
    
    ref_z = np.interp( cumulative , np.cumsum( zh[0] ), 
                       zh[1][0:-1] + 1 )
    #ref_rho = np.interp( cumulative, np.cumsum( rhoh[0] ), 
    #                   rhoh[1][0:-1] + 0.125 ) 
    
    ref_sigma = np.interp( cumulative , np.cumsum( sigmah[0] ), 
                      sigmah[1][0:-1] + 0.125 )
    ref_rho = rho_from_sigma( ref_sigma, ref_z )

    # repackage ref_rho as xr.dataarray
    ref_rho = xr.DataArray( data = ref_rho , coords = \
               {'z_t':ref_z } )
    ref_rho = ref_rho.interp( z_t = ds['z_t'] )
    
    return ref_rho 


