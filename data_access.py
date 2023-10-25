# data_access.py
'''
Module featuring routines that enable easy access to files stored in Taroko, Berimbau, and Kage. 
'''
import xarray as xr; import pandas as pd; 
from datetime import datetime, timedelta; 
import numbers
# Import mapper is missing to run standardize_coords() on dataset.ds_loader(). Add when I sort out conda environment problems. 

# Route to notable SST products.
urls = dict(); 
urls['HadISST'] = 'http://kage.ldeo.columbia.edu:81/home/.datasets/.HadISST/.HadISST_sst.nc/.sst/dods'
urls['ORAs5'] = 'http://kage.ldeo.columbia.edu:81/home/.datasets/.ORAs5/.ingrid-ready/.sst.nc/.sst/dods'
urls['COBE'] = 'http://kage.ldeo.columbia.edu:81/home/.datasets/.COBE-SST/.sst.mon.mean.nc/.sst/dods'
urls['COBE2']= 'http://kage.ldeo.columbia.edu:81/home/.datasets/.COBE-SST2/.sst.mon.mean.nc/.sst/dods'
urls['ERSSTv5'] = 'http://kage.ldeo.columbia.edu:81/home/.datasets/.ERSST/.ERSSTv5.nc/.sst/zlev/removeGRID/dods'
start = {'HadISST':1870,'ORAs5':1958,'COBE':1891,'COBE2':1850,'ERSSTv5':1854}

dt = {'HadISST':timedelta(days = 1), 'ORAs5':timedelta(days=1), 
      'ERSSTv5':timedelta(days= 365.25/12 ),  'COBE':timedelta(days=1), 
      'COBE2': timedelta(days=1) }

def hour_rounder(t):
    # Rounds to nearest hour by adding a timedelta hour if minute >= 30
    return (t.replace(second=0, microsecond=0, minute=0, hour=t.hour)
               +timedelta(hours=t.minute//30))

class dataset:
    # Class for gridded datasets. 
    def __init__( self, name ):
        self.name = name;
        self.url = urls[name];
        self.start = datetime( start[name], 1, 1, 0 ); # hour zero of Jan 1 of given years
        self.dt = dt[name]
        self.ds = self.ds_loader()
        
    def ds_loader( self ):
        temp_ds = xr.load_dataset( self.url, decode_times = False );
        temp_ds = self.time_fixer( temp_ds );
        return temp_ds
    
    def time_fixer( self , xr ):
        try:
            if isinstance( xr['time'].values[0] , numbers.Number):
                xr['time'] = [ hour_rounder( \
                        self.start + self.dt * jj ) for jj in xr['time'].values ]
            else: 
                print( 'Time coordinates are not numeric')
        except:
            print('Coordinate time does not exist.')
        return xr
        