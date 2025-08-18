import xarray as xr; import numpy as np;

# Basic differentiation functions

def deriv_x( xr_obj , ds = None, lat0 = None ):
    # Compute derivative over x (longitude) using meters
    
    if lat0 is None:
        norm = np.cos( xr_obj['lat'] / 180 * np.pi )
    else:
        norm = np.cos( lat0 / 180 * np.pi )

    objx = xr_obj.differentiate( 'lon' ) / ( 110e3 * norm )
    return objx

def deriv_y( xr_obj, ds = None ):
    objy = xr_obj.differentiate( 'lat' ) / 110e3
    return objy

def deriv_z( xr_obj , ds ):
    # Get information about vertical dimension
    #zind = find_dim_index( xr_obj , 'nVertLevels' )
    #dz = np.diff( ds[ 'zMid' ].values , axis = zind )
    dz = ds['zMid'].differentiate( 'nVertLevels' )

    return xr_obj.differentiate( 'nVertLevels' ) / dz


def gradient_vec( xr_obj, ds ):
    # Get gradient of xr_obj, which must be defined on nVertLevels (center of grid cells)
    grad_comps = [ deriv_x( xr_obj, ds ),
                   deriv_y( xr_obj, ds ),
                   deriv_z( xr_obj, ds ) ]
    # Put together as a vector
    for ii in range(3):
        grad_comps[ii] = grad_comps[ ii ].expand_dims( \
                dim = { 'diff_along' : [ ii ] } , axis = 0 )
    # Squeeze into single object
    gradient = xr.concat( grad_comps, dim = 'diff_along' )
    return gradient

# Book keeping functions ( know your dataset )
def find_dim_index( xr_obj, dim_name ):
    try:
        index = xr_obj.dims.index( dim_name )

    except ValueError:
        print( f"'{dim_name}' not found" )

    return index

def uvec( ds ):
    # Return flow components from E3SM ordered as a vector in xarray
    flow_comps = [ ds['velocityZonal'],
                   ds['velocityMeridional'],
                   w_in_center( ds ) ]
    # Add flow direction dimension to each component
    for ii in range(3):
        flow_comps[ii] = flow_comps[ii].expand_dims( \
                      dim = { 'direction' : [ ii ] } , axis = 0 )
    return xr.concat( flow_comps, dim = 'direction' )


def w_in_center( ds ):
    # Given an E3SM dataset, return w at the center
    # of grid cells. 
    
    # Move from one dimension to another (both along same axis)
    orig_dim_name = 'nVertLevelsP1'
    target_dim_name = 'nVertLevels'

    # Get current locations and interpolate to desired
    zloc = ds[ orig_dim_name ].values + 0.5
    w = ds['vertVelocityTop'].interp( { orig_dim_name : zloc } ) # expect nan at bottom

    # Repackage onto target dimension
    target_w = xr.DataArray( data = w.values , dims = ds['velocityZonal'].dims , 
                             coords = ds['velocityZonal'].coords )
    
    return target_w

def strain_rate_tensor( ds ):
    # Input is data from e3sm
    # Get flow components as xarray
    u_i = uvec( ds )

    # Define methods for differentiation
    diff_funcs = [ deriv_x, deriv_y, deriv_z ]
    grads = [ ]
    
    for jj in range(3):
        # Cycle through dimensions of differentiation
        differ = diff_funcs[ jj ]
        grads.append( differ( u_i , ds ).expand_dims( \
                        dim = { 'diff_along' : [ jj ] }, axis = 1 ) )

    # Concatenate along j
    djui = xr.concat( grads, dim = 'diff_along' ); 
    # Transpose i and j
    tdims = list( djui.dims );
    # Use indices to reorder old dimensions
    complem = [ kk for kk in range( 2, len( djui.shape ) ) ]
    nudims = tuple(  tdims[kk]  for kk in ( [1,0] + complem )  )
    # Transpose onto the reordered dimensions
    diuj = djui.transpose( *nudims )
    Sij = ( djui + diuj.values ) / 2
    return Sij





