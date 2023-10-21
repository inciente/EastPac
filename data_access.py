# data_acces.py
'''
Module featuring routines that enable easy access to files stored in Taroko, Berimbau, and Kage. 
'''
import xarray as xr; import pandas as pd; 
from datetime import datetime, timedelta; 

# Route to notable SST products.
urls = dict(); 
urls['HadISST'] = 'http://kage.ldeo.columbia.edu:81/home/.datasets/.HadISST/.HadISST_sst.nc/.sst/dods'
urls['ORAs5'] = 'http://kage.ldeo.columbia.edu:81/home/.datasets/.ORAs5/.ingrid-ready/.sst.nc/.sst/dods'
urls['COBE'] = 'http://kage.ldeo.columbia.edu:81/home/.datasets/.COBE-SST/.sst.mon.mean.nc/.sst/dods'
urls['COBE2']= 'http://kage.ldeo.columbia.edu:81/home/.datasets/.COBE-SST2/.sst.mon.mean.nc/.sst/dods'
urls['ERSSTv5'] = 'http://kage.ldeo.columbia.edu:81/home/.datasets/.ERSST/.ERSSTv5.nc/.sst/zlev/removeGRID/dods'
start = {'HadISST':1870,'ORAs5':1958,'COBE':1891,'COBE2':1850,'ERSSTv5':1854}

class dataset:
    # Class for gridded datasets. 
    def __init__( self, name ):
        self.name = name;
        self.url = urls[name];
        self.start = datetime( start[name], 1, 1, 0 ); # hour zero of Jan 1 of given years
        
    
    @property
    def ds(self):
        return self._ds

    @ds.setter
    def ds( self, opts = None ):
        temp_ds = xr.load_dataset( self.name , decode_times = False );
        self._ds = time_maker( temp_ds );






