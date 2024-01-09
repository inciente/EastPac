import xarray as xr; from scipy import linalg
import numpy as np;
import pandas as pd; from datetime import datetime, timedelta;

def seasonal_matrix( timevec, n_terms , leap_years = False ):
    '''
    Return a matrix with mean, linear trend, and n_terms pairs of cosines and cosines
    meant to be used with least squares fits. Frequencies predetermined to be a seasonal cycle.
    '''
    year = 365
    if leap_years:
        year = 356.25

    N = len( timevec ); 
    # t0 must be middle point in order for coef0 to be the mean
    t0 = datetime( 1900, 1, 1 )
    xaxis = np.array([ (pd.to_datetime(jj) - t0 ).total_seconds() \
            for jj in timevec.values ])
    year_freq = 2 * np.pi / (3600*24*year)

    # Create empty matrix
    fit_mat = np.empty( (N, 2*n_terms + 2 ) ); # mean, trend, cycles
    fit_mat[:,0] = np.ones( (N, ) ); 
    fit_mat[:,1] = xaxis; # linear fit

    # Put sines and cosines into matrix
    for jj in range( 1, n_terms + 1 ):
        cols = [jj*2, jj*2 + 1];
        freq = year_freq * jj 
        fit_mat[:,cols[0]] = np.sin( xaxis * freq )
        fit_mat[:,cols[1]] = np.cos( xaxis * freq )
    return fit_mat 



def seasonal_cycle( data, n_terms, leap_years = True ):
    ''' 
    Use least squares to fit a trend and seasonal cycle to data. 
    '''
    if len( data.values.flatten().shape ) > 1:
        print('I can only process one-dimensional data')
        return None

    # Create sine cosine fit matrix
    fit_mat = seasonal_matrix( data['time'], n_terms, leap_years )

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


def lagged_correlation(x, y, max_lag=10):
    """
    Compute lagged correlation between two xr.DataArrays.

    Parameters:
    - x, y: xr.DataArray
        Input arrays with the same dimensions.
    - max_lag: int, optional
        Maximum lag to consider along the time dimension.

    Returns:
    - lagged_corr: xr.DataArray
        Lagged correlation along the time dimension.
    """
    # Check if input arrays have the same dimensions
    if x.dims != y.dims:
        raise ValueError("Input arrays must have the same dimensions.")
    # Initialize an array to store lagged correlations
    lagged_corr_values = np.zeros((2 * max_lag + 1,) + x.shape[1:], dtype=np.float64)

    # Iterate over lags
    for lag in range(-max_lag, max_lag + 1):
        # Apply lag along the time dimension
        x_lagged = np.roll(x.values, lag, axis=0)
        y_trimmed = y.values[max(0, -lag):, ...]
        x_trimmed = x_lagged[max(0, lag):, ...]

        # Compute correlation along all dimensions
        corr_values = np.corrcoef(x_trimmed.flatten(), y_trimmed.flatten())[0, 1]
        lagged_corr_values[lag + max_lag] = corr_values

    # Create a DataArray for lagged correlations
    lagged_corr_dims = ["lag"] + list(x.dims[1:])
    lag_values = np.arange(-max_lag, max_lag + 1, dtype=np.int64)
    lagged_corr = xr.DataArray(lagged_corr_values, dims=lagged_corr_dims, coords={"lag": lag_values})

    return lagged_corr


def get_trend( xr_obj, for_dims ):
    # get temporal trend (per decade) in xr_obj
    # over_dims is a list indicating all non-temporal dimensions
    xr_obj = xr_obj.stack( all_points = over_dims )
    timevec = np.array( [ (kk - datetime(1980,1,1) ).total_seconds() \
                for kk in pd.to_datetime( xr_obj['time'].values ) ] )
    # replace nans with zeros
    mask = np.isnan( xr_obj.values );# mask = np.isnan( vals ); 
    #vals[ mask ] = 0; # because polyfit crashes with nans
    coefs = np.polyfit( timevec, xr_obj.values, 1 ) 
    # put in mask for nans
    coefs[ : , mask[0,:] ] = np.nan
    # save into
    nu_dat = coefs[0] * (10*365*24*3600 ); # ttrend per decade
    # save as xr
    coords = xr_obj['all_points'].coords;
    nu_obj = xr.DataArray( data = nu_dat , coords = coords ).unstack('all_points' )
    return nu_obj  




