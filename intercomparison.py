#navigate.py

from abc import ABC, abstractmethod
import sys, os, re;
import xarray as xr;
import pandas as pd; 

def do_nothing( data_in ):
    # Is useful as placeholder in some more complex functions
    return data_in

class intercomparison:
    '''
    Routines that will help access large model output within complex folder structures.
    - - - To set up a comparison, you'll need a dir_table, which includes information about the model runs and respective directories. 
    - - - Here's an example of what a dir_table (pandas dataframe) may look like: 
    ----- |  config  | ensnum |     folder      | 
    ----- | control  |   01   |  f01_g12.H.c1   |
    ----- | fluxadj  |   01   |  f01_g12.H_fa1  |
    '''
    def __init__( self , dir_table ):
        # Properties including paths, config, etc.
        self.dir_table = dir_table; 
        self.models = self.make_models();

    def make_models( self ):
        models = []; 
        for jj in range( len( dir_table.index ) ):
            # Create instance of model_run for each row
            models.append( model_run( dir_table[jj] ) )
    
    

class model_run(ABC):
    '''
    Subclass that helps compose a model_comparison object. 
    '''
    def __init__(self,  dir_row ):
        self.config = dir_row['config'];
        self.path = dir_row['path'];
        self.ensnum = dir_row['ensnum'];
        #self.components = folders_within(); # subdirectories inside main folder
        

    def folders_within( self ):
        # Get a list of all subdirectories within this model run folder
        dl_clean = [];
        for dir_name in os.listdir( self.path ):
            # Check what items within are directories
            if os.path.isdir( os.path.join( path, dir_name )):
                dl_clean.append( dir_name );
        return dl_clean         

    def files_rule( self, subdir, rule):
        # Get all files within subdir, and return a list of all files within whose naming follows a given rule. 
        all_files = os.listdir( self.path + subdir );
        is_file_good = [ rule( item ) for item in all_files ]; # rule returns true or false 
        # Now save only good files and add source path so we can load them easily
        good_files = [ self.path + subdir + '/' + all_files for all_files, is_file_good \
                     in zip( all_files, is_file_good ) if is_file_good ]; 
        good_files = sorted( good_files ); 
        return good_files 

    @abstractmethod
    def load_batch( self ):
        # This will be a very important function. Allow access to a subset of the files in file_list. 
        # Enable iteration in cycle that includes all our actual operations and analysis on model output. 
        pass


    @abstractmethod
    def prepare_xr( self ):
        # Here is where to code in specific requirements for coordinate assignment/standardization, etc. 
        pass


 
class CESM_12(model_run):
    '''
    Originally made to analyze Boniface's runs of GHG emissions under flux adjustment. 
    '''

    def prepare_xr( self, xr_obj ):
        return xr_obj 

    def load_batch( self , filelist, cut_func=do_nothing ):
        # Take a list of files, load all of them, and concatenate over time
        # To save storage, cut files in space or subselecting variables using cut_func
        ind_files = []; 
        filelist = sorted( filelist ) ; # sort it just in case 

        for jj in range( len( filelist ) ):
            filename = filelist[jj];
            datafile = xr.load_dataset( filename , decode_times = False )
            # Now prepare xr to have desired format and cut it spatially /select variables as desired
            datafile = self.prepare_xr( datafile );
            datafile = cut_func( datafile );
            ind_files.append( datafile ); 

        del datafile # keep it clean
        batch_file = xr.concat( ind_files , dim = 'time' ); 
        return batch_file 






