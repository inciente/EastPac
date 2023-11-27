#navigate.py

from abc import ABC, abstractmethod
import sys, os, re, time, warnings;
import xarray as xr;
import pandas as pd; 

def do_nothing( data_in ):
    # Is useful as placeholder in some more complex functions
    return data_in

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
        self.storage = self.prepare_storage();

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

    def prepare_storage( self ):
        # Create dict with entries for each one of the configurations in dir_table.
        storage = dict(); 
        configs = list( set( self.dir_table['config'] ) ) 
        for conf in configs:
            storage[conf] = dict()
        return storage


    def distribute_task( self, to_storage = False ):
        # Cycle through models, extract data, perform task, and save if necessary 
        pass



class model_run(ABC):
    '''
    Subclass that helps compose a model_comparison object. 
    '''
    chunks = dict()
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
        #subdir = self.comp['subdir']; # -- don't include comp here
        #rule = self.comp['file_rule']; 
        # Get all files within subdir, and return a list of all files within whose naming follows a given rule. 
        all_files = os.listdir( self.path + subdir );
        is_file_good = [ rule( item ) for item in all_files ]; # rule returns true or false 
        # Now save only good files and add source path so we can load them easily
        good_files = [ self.path + subdir + '/' + all_files for all_files, is_file_good \
                     in zip( all_files, is_file_good ) if is_file_good ]; 
        good_files = sorted( good_files ); 
        return good_files 

    def mf_loader( self, filelist, cutter = do_nothing ):
        # General function used to load and concatenate files. dropper allows to drop useless or confounding variables
        # Writing it here as a placeholder for whenever I figure out the best way to load large numbers of files
        ind_files = [];
        print( filelist )
        filelist = sorted( filelist ); # just in case
        batch_file = xr.open_mfdataset( filelist, chunks = {'time':1} , parallel = False, combine='by_coords');
        batch_file = self.prepare_xr( batch_file )
        batch_file = cutter( batch_file );
        return batch_file

    #@abstractmethod
    #def load_batch( self ):
        # This will be a very important function. Allow access to a subset of the files in file_list. 
        # Enable iteration in cycle that includes all our actual operations and analysis on model output. 
    #    pass


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
    def prepare_xr( self, xr_obj ):
        '''
        Focus on changing coordinate names and sorting by ascending order
        '''
        xr_obj['z_w_top'] = xr_obj['z_w_top']/100; 
        xr_obj['z_t'] = xr_obj['z_t'] / 100; # cm to meters
        xr_obj["nlon"] = xr_obj["ULONG"].isel(nlat=0);
        xr_obj["nlat"] = xr_obj["ULAT"].isel(nlon=0);
        xr_obj = xr_obj.sortby( 'nlon' ); 
        # Interoperability with timedelta
        xr_obj['time'] = xr_obj.indexes['time'].to_datetimeindex()
        return xr_obj 

    #def load_batch( self , filelist, cut_func = do_nothing ):
    #    # Take a list of files, load all of them, and concatenate over time
    #    # To save storage, cut files in space or subselecting variables using cut_func
    #    ind_files = []; 
    #    filelist = sorted( filelist ) ; # sort it just in case 
    #    batch_file = xr.open_mfdataset( filelist , chunks = {'time':12}, parallel = False,
    #                          combine = 'by_coords' );
        #popem = True
        #if popem : 
        #    dum = xr.open_dataset( filelist[0] )
        #    batch_file.update( dum[[ 'ULONG','ULAT','TLONG','TLAT' ]] )
        #    dum.close()
        # Fix coordinate situation
    #    batch_file = self.prepare_xr( batch_file );
    #    batch_file = cut_func( batch_file ); # get from intercomparison.scope
    #    return batch_file 

class CAM( model_run ):
    '''
    Set of functions needed to access data within the atmospheric output of CESM 1.2, 
    coming from the Community Atmosphere Model (CAM).
    '''
    chunks = {'time':1};
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
