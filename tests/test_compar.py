'''
Unit and integration tests for the code in intercomparison.py
'''

import unittest, os, sys;
import xarray as xr; import pandas as pd; 
import numpy as np

# Route from tests folder to EastPac modules
test_folder = os.getcwd()
EP_folder = os.path.split( test_folder )[0];
# Import intercomparison module for testing
sys.path.append( EP_folder )

import intercomparison
#target = __import__( 'intercomparison.py' );

'''
Create files and a pd.DataFrame to execute tests.
'''

latitudes = np.array([-2,-1,0,1,2]); 
times = np.arange( 2, 100, 5) + 5e5; # four files, days after 0000-01-01

def create_random_xr( lat, tval ):
    coords = { 'time' : np.array( [tval] ), 'lat' : lat }
    data_vals = np.expand_dims( np.random.rand( len( lat ) ), 0 ); 
    rand_data = xr.DataArray( data = data_vals, coords = coords, name = 'dummy' ); 
    return rand_data 

for jj in range( len( times ) ):
    # Create one file for each time - for model configuration 1
    data_bit = create_random_xr( latitudes, times[jj] ); 
    # Name of file under which this will be stored
    file_bit = 'config01/dat/data_test_' + str(jj).zfill(2) + '.nc'; 
    data_bit.to_netcdf( file_bit ); # create individual netcdf file
    del data_bit

    # Now repeat the operation for model configuration 2 
    data_bit = create_random_xr( latitudes, times[jj] );
    file_bit = 'config02/dat/data_test_' + str(jj).zfill(2) + '.nc';
    data_bit.to_netcdf( file_bit ); # create another netcdf file
    del data_bit

# Create dataframe that will allow to load these files into model_run subclasses
my_runs = pd.DataFrame.from_dict( {'config':['standard','standard'], 'path':['config01/','config02/'], 
        'ensnum':[1,2] } );

# Create intercomparison and instantiate models
comp = intercomparison.intercomparison( my_runs, [intercomparison.DUMMY]*len(my_runs) )
# Give comparison object a scope 
scope = dict();
scope['load_from'] = 'dat'
scope['file_rule'] = lambda text: True;#lambda fln : True if '5-01.basic.nc' in fln else False
scope['cutter'] = lambda xr_obj: xr_obj; # return xarray without changing
comp.scope = scope;
print('Created intercomparison, instantiated models, and assigned a scope')

ds = comp.extract_data(0)
print('Successfully ran extract data') 
print( ds )

# All tests for intercomparison classes
class TestComparison( unittest.TestCase ):

    

    def test_model_creation( self ):
        self.assertGreater( 4 , 1);
        
    def find_files( self ):
        self.assertEqual( 2, 2)                

        def my_rule( filename ):
            return ( 'data_test' in filename )
        


