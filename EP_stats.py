import xarray as xr; from scipy import linalg
import numpy as np;
import pandas as pd; from datetime import datetime, timedelta;



def seasonal_cycle( data, n_terms ):
    ''' 
    Use least squares to fit a trend and seasonal cycle to data. 
    '''
    if len( data.values.flatten().shape ) > 1:
        print('I can only process one-dimensional data')
        return None

    # Specify parameters needed to create sines, cosines, etc.
    N = len(data); # length of record
    t0 = datetime(1900,1,1)
    xaxis = np.array([ (pd.to_datetime( jj )  - t0 ).total_seconds() \
            for jj in data.time.values ]); # for trend
    year_freq = 2 * np.pi / (3600*24*365.25);

    # Create empty matrix for base elements of fit
    fit_mat = np.empty( ( N, 2*n_terms + 2 ) ); # mean, trend, and cycles
    fit_mat[:,0] = np.ones( ( N , ) );
    fit_mat[:,1] = xaxis

    # Put sines and cosines into matrix
    for jj in range(1,n_terms+1):
        cols = [jj*2, jj*2+1]
        freq = year_freq * jj
        fit_mat[:, cols[0]] = np.sin( xaxis * freq );
        fit_mat[:, cols[1]] = np.cos( xaxis * freq );

    # Obtain coefficients for all terms
    coefs, res = solve_least_squares( fit_mat, data )
    # Create new time series to match original data
    nu_data = np.matmul( fit_mat, coefs )
    nu_data = xr.DataArray( data = nu_data, coords = data.coords );

    return coefs, nu_data
        
        

def solve_least_squares( fit_mat, data ):
    '''
    Simple linear algebra solution to fit linear coefficients to columns in fit_mat, such that 
    least squares misfit is minimized. 
    '''
    coefs, res, rnk, s = linalg.lstsq( fit_mat, data  )  
    return coefs, res


