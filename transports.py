import xarray as xr; import numpy as np; 
import pandas as pd; import sys, os, time;
import shapely as shp

''' 
Suite of tools to calculate transports across transects, into and out of boxes. 

'''

class transect:
    def __init__( self, x, y ):
        pass

def normal_vec( line ):
    # Get vector representing normal orientation to line
    pass

def flux_across( vec, line ):
    # Project vectors onto the normal direction to line
    # line should be instance of transect
    norm_dir = normal_vec( line );
    locs = [1,1,1]; # set of points along line
    flux = 0;
    for jj in range( locs ):
        vec_here = vec.sel( at_loc[jj] );
        flux_here = np.dot( vec_here, norm_dir ); # normal component
        flux += flux_here * dl
    pass

class ocean_boxes:
    def __init__( self, polys ):
        self.polys = pol_div; # instance of shp.geometry.multipolygon

    


