# intercomparison.py

from datetime import timedelta, datetime
from abc import ABC, abstractmethod
try:
    from kerchunk.hdf import SingleHdf5ToZarr
    from kerchunk.netCDF3 import NetCDF3ToZarr
    from kerchunk.combine import MultiZarrToZarr
    import ujson
except:
    print('Could not import kerchunk. Implementations in intercomparison will fail')
#from datatree import DataTree

import sys, os, re, time, warnings;
import xarray as xr;
import pandas as pd; 
import fsspec;

def do_nothing( data_in ):
    # Is useful as placeholder in some more complex functions
    return data_in

def days_since_to_datetime( t0, timevec, no_leap = True ):
    # Input: t0 -> datetime, timevec -> array of numbers
    # Output: time -> array of datetimes
    ylen = 365.25
    if no_leap:
        ylen = 365;
    years = [ np.floor( tval / ylen ) + t0.year for tval in timevec ] ;
    ydays = [ np.mod( tval, ylen ) for tval in timevec ] ;
    time = np.array( [ datetime( years[jj] ) + timedelta( days = ydays[jj] ) \
                for jj in range( len( timevec ) ) ] )
    return time

def write_json( fs_from, fs_to, filepath, write_as ):
    '''
    Function to take a single netcdf file and make json file using kerchunk 
    Input: 
    --- filepath --- str leading to a netcdf file
    --- fs_from --- fsspec.filesystem from which to extract data
    --- fs_to --- fsspec.filesystem to write json files with
    Output:
    --- Writes json file that allows xr to access to netcdf as if it were a zarr    
    '''

    so = dict( mode = 'rb', anon = True, default_fil_cache = False, 
            default_cache_type = 'first' ); # args to fs.open()
    #with fs_from.open( filepath , **so ) as infile: 
    #h5chunks = SingleHdf5ToZarr( infile, filepath, inline_threshold = 350 ); 
    h5chunks = NetCDF3ToZarr( filepath , inline_threshold = 300 ); 
    # inline_threholds is directly proportional to size of json, inverse to time it takes to open data later
    with fs_to.open( write_as, 'wb' ) as writer:
        writer.write( ujson.dumps( h5chunks.translate()).encode() );

def json_combiner( fs_to, files2combine, config, save_as ):
    # combine all individual jsons to create a wormhole to whole dataset
    mzz = MultiZarrToZarr( files2combine , remote_protocol = 'file',
                      concat_dims = ['time'],  
                      identical_dims = config['id_dims'] )
    d = mzz.translate()
    print('Combining ' + str( len(files2combine) ) + ' different files!')
    with fs_to.open( save_as , 'wb' ) as f:
        f.write( ujson.dumps( d ).encode())
    
def json_to_xr( fln ):
    #fln = get_summary_fln( model ); # get filepath
    # Load dataset from json as xarray 
    fs = fsspec.filesystem('reference', fo= fln, target_options={'anon':True}, remote_options={'anon':True})
    m = fs.get_mapper('')
    ds = xr.open_dataset(m, engine='zarr', backend_kwargs={'consolidated':False}, 
            chunks = {} )
    #ds = model.prepare_xr( ds ) # remnant from implementation with model_run instances
    return ds 



class ScopeDescriptor:
    '''
    Use this class (dictionary-like value) to store all information needed to subset, 
    modify, acces, or generally operate with entities handler by model_run instances.
    This is a sort of interface for intercomparison to know what to say to model_run.
    '''

    def __get__( self, instance, owner ): 
        return instance._scope

    def __set__( self, instance, value ):
        # Make all sanity checks needed to ensure that scope will work seamlessly
        if not isinstance( value, dict ):
            raise ValueError('Scope must be a dictionary-like object' )
        # Add other sanity checks to make sure that scope is well-defined
        if not all( key in value.keys() for key in ['load_from','cutter','file_rule'] ):
            incomplete_scope(); # tell user about missing entries
        instance._scope = value
    
    def incomplete_scope( ):
        warning_message = 'Your scope does not include entries for: load_from, cutter, and file_rule. Do you want to continue? (y/n):'
        user_input = input( warning_message ).lower()
        if user_input != 'y':
            raise UserWarning('Execution stopped by user')
        else:
            warnings.warn('User chose to proceed with the specified scope', UserWarning )



class intercomparison:
    '''
    Routines that will help access large model output within complex folder structures.
    - - - To set up a comparison, you'll need a dir_table, which includes information about the model runs and respective directories. 
    - - - Here's an example of what a dir_table (pandas dataframe) may look like: 
    ----- |  config  | ensnum |     folder      | 
    ----- | control  |   01   |  f01_g12.H.c1   |
    ----- | fluxadj  |   01   |  f01_g12.H_fa1  |
    '''
    
    def __init__( self , dir_table, model_types ):
        # Properties including paths, config, etc.
        self.dir_table = dir_table ; 
        # iterable with subclass of model_run should be used for each model
        self.model_types = model_types; # could eventually replace with factory managed by scope
        # create instance for each row in the model table
        self.models = [ model_types[jj]( dir_table.loc[jj] ) for jj in range( len( dir_table ) ) ]; 

    # Use the descriptor for the scope attribute
    scope = ScopeDescriptor()

    def extract_data( self, m_index ):
        # Basic routine of operations needed to prepare, load, and cut data as required by scope 
        model = self.models[ m_index ] 
        # create a list of files to load
        files2load = model.files_rule( self.scope['load_from'], self.scope['file_rule'] )
        # now use load_batch in the model to load them
        start_time = time.time()
        data = model.mf_loader( files2load, self.scope['cutter'] ); 
        print("Loading batch took --- %s seconds ---" % (time.time() - start_time))
        # data = data.sortby('time');
        return data 

    def model_to_json( self, m_index, override_exists = False ):
        '''
        Takes in data from scope and model to write json summary files of multiple .nc

        '''
        model = self.models[ m_index ]; 
        files2write = model.files_rule( self.scope['load_from'] , self.scope['file_rule'] ); 
        # Cycle through individual files and write each a json
        files2combine = []
        for jj in range( len( files2write ) ):
            get_from = files2write[jj];
            save_as = self.scope['save_as']( model, get_from );
            files2combine.append( save_as )
            if os.path.isfile( save_as ) and override_exists is not True:
                continue
            else:
                # access netcdf, open, and write json 
                write_json( self.scope['fs_from'], self.scope['fs_to'], get_from, save_as )
        return files2combine


    def prepare_storage( self ):
        # Create dict with entries for each one of the configurations in dir_table.
        storage = dict(); 
        configs = list( set( self.dir_table['config'] ) ) 
        for conf in configs:
            storage[conf] = dict()
        return storage

    def stack_config( self, conf_str ):
        # Returns a list with model instances under a given configuration
        sub_table = self.dir_table.loc[ self.dir_table['config'] == conf_str ]
        # Create list with model instances in sub_table
        models = [ self.models[kk] for kk in sub_table.index ] 
        return models

    def create_comp_xr( self ):
        # Create dataset with dims necessary to store output from all models,
        # ordered by config and ensnum
        pass 

class model_run(ABC):
    '''
    Subclass that helps compose a model_comparison object. 
    '''
    # Identifying features of model instances
    chunks = dict()
    name = str()
    mzz_config = dict()

    def __init__(self,  dir_row ):
        self.config = dir_row['config'];
        self.path = dir_row['path'];
        self.ensnum = dir_row['ensnum'];
        #self.components = folders_within(); # subdirectories inside main folder        

    def view_subdirs( self ):
        # Get a list of all subdirectories within this model run folder
        dl_clean = [];
        for dir_name in os.listdir( self.path ):
            # Check what items within are directories
            if os.path.isdir( os.path.join( self.path, dir_name )):
                dl_clean.append( dir_name );
        return dl_clean         

    def files_rule( self, subdir, rule ):
        # Get all files within subdir, and return a list of all files within whose naming follows a given rule. 
        # Tool needed to extract files
        fs = fsspec.filesystem('', anon=True);
        # Where we're going to look
        folder = self.path + subdir + '/';
        # Find files that look like the rule
        good_files = sorted( fs.glob( folder + rule ) )
        return good_files 

    def mf_loader( self, filelist, cutter = do_nothing ):
        # General function used to load and concatenate files. dropper allows to drop useless or confounding variables
        # Writing it here as a placeholder for whenever I figure out the best way to load large numbers of files
        batch_file = xr.open_mfdataset( filelist, chunks = {'time':1, 'nlat':1} , parallel = True, combine='by_coords');
        batch_file = self.prepare_xr( batch_file )
        batch_file = cutter( batch_file );
        return batch_file

    @abstractmethod
    def prepare_xr( self ):
        # Here is where to code in specific requirements for coordinate assignment/standardization, etc. 
        pass



 
class POP2(model_run):
    '''
    Oriiginally made to analyze runs of GHG emissions under flux adjustment. 
    Current version can only handle ocean component. 
    '''
    chunks = {'time':1, 'nlon':6}
    name = 'POP2'
    mzz_config = { 'id_dims':['nlon','nlat','z_t','z_w_top'] }

    def prepare_xr( self, xr_obj ):
        '''
        Focus on changing coordinate names and sorting by ascending order
        '''
        for vvar in ['z_t','z_w_top','z_w_bot']:
            try:
                xr_obj[ vvar ] = xr_obj[ vvar ] / 100; # scaled to m
            except:
                pass
        xr_obj["nlon"] = xr_obj["ULONG"].isel(nlat=0).isel( time = 0 );
        xr_obj["nlat"] = xr_obj["ULAT"].isel(nlon=0).isel( time = 0 );
        #elif isinstance( xr_obj, xr.DataArray ):
        xr_obj = xr_obj.sortby( 'nlon' ); 
        xr_obj = xr_obj.rename( {'nlon':'lon', 'nlat':'lat'} )
        # Interoperability with timedelta
        # substitute for time fixer days since 0000-01-01
        xr_obj['time'] = xr_obj.indexes['time'].to_datetimeindex()
        return xr_obj 


class CAM( model_run ):
    '''
    Set of functions needed to access data within the atmospheric output of CESM 1.2, 
    coming from the Community Atmosphere Model (CAM).
    '''
    chunks = {'time':1};
    name = 'CAM'
    mzz_config = {'id_dims':['lat','lon','lev']}

    def prepare_xr( self, xr_obj ): 
        '''
        Get values on most workable format
        '''
        xr_obj['time'] = xr_obj.indexes['time'].to_datetimeindex();
        # As far as I can tell, everything else has a standard name (lon, lat, time, lev) and value
        return xr_obj 

    #def load_batch( self, filelist, cut_func = do_nothing ):
    #    pass



class DUMMY( model_run ):
    '''
    Not a real model, just a dummy version designed to test operation of intercomparison
    '''
    chunks = {'time':1};
    def prepare_xr( self, xr_obj ):
        return xr_obj 
