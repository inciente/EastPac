import xarray as xr; from scipy import linalg
import numpy as np;
import pandas as pd; from datetime import datetime, timedelta;
from scipy import signal

def numeric_time( xr_time ):
    # Input is xr.DataArray representing time coordinate
    # Output must be same data in numeric format, with units of seconds.
    t0 = datetime( 1900, 1, 1);
    timevec = np.array( [ ( pd.to_datetime( tt ) - t0 ).total_seconds() \
                         for tt in xr_time.values ] );
    return timevec

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
    xaxis = numeric_time( timevec )
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

def remove_seasonal( xr_obj , leap_years = False ):
    # Get detrended anomalies of xr_obj without seasonal cycle
    ncomp = 4; # harmonic components for seasonal cycle
    coefs, recreated = seasonal_cycle( xr_obj, ncomp, 
             leap_years = leap_years )
    return xr_obj - recreated


def fancy_lagged_correlation( x, y , max_lag = None, dim = 'time'):
    # Takes in two xr time series and computes lagged correlation with 1-sized steps along dim
    if max_lag is None:
        max_lag = int( len( x[ dim ] ) * 0.75 )
    corr_vals = [ xr.corr( x, y, dim = dim ) ] \
                 + [ xr.corr( x[n:], y[:-n] ) for n in range(1, max_lag + 1 ) ]
    corr_inverse = [ xr.corr( y[-m:], x[:m] ) for m in range(-max_lag, 0 ) ]
    all_corr = corr_inverse + corr_vals
    lags = [ jj for jj in range( -max_lag, max_lag + 1 ) ] 
    all_corr = xr.DataArray( data = all_corr, dims = ('lag'), coords = {'lag':lags} )
    return all_corr

def lagged_correlation(x, y, max_lag = None):
    # Computed lagged correlations between two xr_objs or 1D arrays
    normalize = lambda data : ( data - np.mean( data ) ) \
                              / np.std( data )
    # This current version computes normal correlation, does not take care of lagging
    c = np.correlate( normalize( x ) / len( x ) , 
                    normalize( y ) , 'full' );
    return c


def get_trend( xr_obj, for_dims ):
    # get temporal trend (per decade) in xr_obj
    # over_dims is a list indicating all non-temporal dimensions
    is_1D = len( for_dims ) == 0;
    if not is_1D:
        # unless object is already 1D
    	xr_obj = xr_obj.stack( all_points = for_dims )
    timevec = numeric_time( xr_obj['time'] )

    # replace nans with zeros
    mask = np.isnan( xr_obj.values );# mask = np.isnan( vals ); 
    vals = xr_obj.values;
    vals[ mask ] = 0; # because polyfit crashes with nans
    coefs = np.polyfit( timevec, vals , 1 ) 
    
    if not is_1D:
        # put in mask for nans
        coefs[ : , mask[0,:] ] = np.nan
        coords = xr_obj['all_points'].coords
    else:
        coords = None
    # save into
    nu_dat = coefs[0] * ( 10*365*24*3600 ); # ttrend per decade
    # save as xr
    nu_obj = xr.DataArray( data = nu_dat , coords = coords )
    if not is_1D:
        nu_obj = nu_obj.unstack('all_points' )
    return nu_obj  

classic_months = [ [12,1,2], [3,4,5], [6,7,8], [9,10,11] ]
def separate_seasons( xr_obj , months = classic_months ):
    # Group all data together by sets of months indicated in months (these define the seasons)
    by_months = xr_obj.groupby('time.month')
    by_season = [];
    # Now iterate through the seasons
    for jj in range( len( months ) ):
        months_in_season = months[jj]; 
        # Get time indices for data corresponding for this season:
        month_inds = np.concatenate( [ by_months.groups[ mm ] for mm in months_in_season ] ); 
        # Save corresponding data
        by_season.append( xr_obj.isel( time = month_inds ).sortby('time') )
    # Now concatenate onto a single xr object
    by_season = xr.concat( by_season , dim = 'season' )
    return by_season 


















