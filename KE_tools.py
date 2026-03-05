import sys, gsw 
import matplotlib.pyplot as plt; 
import numpy as np; 
import pandas as pd; import xarray as xr; 

# Modules from this package
import vectors

''' 
Functions to compute the energy budget of the upper Pacific Ocean
E3SM branch

'''

def buoyancy_from_rho( rho ):
    return 9.81 * ( rho - 1024 ) / 1024

def flip_profile( xr_obj ):
    # Take a 1D xr object, and:
    # ---  turn its coordinate into its variable
    # ---  turn its variable into coordinate
    coord = list( xr_obj.coords.keys() )[0]
    nu_obj = xr.DataArray( data = xr_obj[ coord ] ,
                          dims = ('variable'), coords = {'variable':xr_obj.values})
    return nu_obj

def get_sigma_ref( ds ):
    # For use with E3SM
    
    # Calculate cell volumes and corresponding densities
    cell_volumes = ( ds['area'] ) * ds['layerThickness']
    all_sigmas = ds['potentialDensity']
    # Flatten data for ease
    flat_sigma = all_sigmas.values.flatten()
    flat_vols = ( cell_volumes ).values.flatten()
    flat_z = ( ds['zMid'] ).values.flatten()

    valid_rows = ( ~ np.isnan( flat_vols ) ) * ( ~np.isnan( flat_sigma ) )
    
    # Get cdf for space available in the ocean space represented by ds
    dz = 5; min_z = np.min( flat_z[ valid_rows ] ) - 2.5;
    n, bins, _ = plt.hist( flat_z[ valid_rows ] , weights = flat_vols[ valid_rows ], 
                                    bins = np.arange( min_z , -2.4, dz ), density = True, cumulative = True )
    #print(n)
    # Organize into a small xr obj
    space = xr.DataArray( data = bins[:-1] + dz / 2 , dims = ('cdf'), 
                                   coords = { 'cdf' : n } ).groupby('cdf').mean()
    
    
    # Get cdf for density of water masses in ds
    drho = 0.1; min_rho = 1021;#np.min( flat_sigma[ valid_rows ] ) + 0.05;
    n, bins, _ = plt.hist( flat_sigma[ valid_rows ] , weights = flat_vols[ valid_rows ], 
                                    bins = np.arange( min_rho, 1035.01 , drho ), 
                                    density = True, cumulative = True )
    #print(n)
    mass = xr.DataArray( data = bins[:-1] + drho / 2 , dims = ('cdf') , 
                                  coords = { 'cdf' : 1 - n } ).groupby( 'cdf' ).mean()
    
    # Now prepare to get corresponding values of z and density
    cdf_vals = np.arange(0.01, 0.9991, 0.01); 
    # now create sigma_ref
    dens_vals = mass.interp( cdf = cdf_vals ).values; 
    z_levels =  space.interp( cdf = cdf_vals ).values
    sigma_ref = xr.DataArray( data = np.concatenate( [ dens_vals, [ dens_vals[-1] ] ] ) , dims = ('z'), 
                             coords = { 'z' : np.concatenate( [ z_levels, [0] ] )} )
    plt.close()
    # get rid of nan values
    sigma_ref = sigma_ref[ ~np.isnan( sigma_ref['z'] ) ]
    # add extra level at surface with lightest density
    single_z = ds['zMid'].sel( {'lon':180, 'lat':0}, method = 'nearest' )
    sigma_ref = sigma_ref.interp( z = single_z , method = 'linear', 
                                 kwargs = {'fill_value' : 'extrapolate' } )
    sigma_ref = sigma_ref.rename( {'nVertLevels':'z'} ); 
    return sigma_ref 



class Tailleux:

    '''
Methods designed to compute most variables needed under the Tailleux formalism. 

Input:
    ds    -    snapshot of POP2 dataset (no time dependence) 
    rho_ref    reference density profile (may be computed from thermo.reference_density() )

All quantities will be computed for all spatial domain in ds, so subset ds before instantiating this class.
    '''
    def __init__( self, ds, sigma_ref ):

        self.ds = ds;
        self.sigma_ref = sigma_ref;
        # save rho values, since they will be used again and again
        self.sigma = ( self.ds['potentialDensity'] ).persist()
        self.nan_mask = ~ np.isnan( self.sigma ); # useful to have
        # compute zr, since it will also be needed repeatedly
        # compressibility not needed since E3SM analysis is based on sigma and not rho
        #self.compressibility = 1024e4 * gsw.kappa( self.ds['SALT'], self.ds['TEMP'], self.ds['z_t'] ).persist()
        self.zr = self.get_zr().persist()

    @property
    def ds(self):
        return self._ds

    @ds.setter
    def ds(self, data):
        if 'Time' in list( data.dims ):
            raise Exception('Cannot instantiate class. Dataset has time dependence.')
        
        self._ds = data

    @property
    def sigma_ref( self ):
        return self._sigma_ref

    @sigma_ref.setter
    def sigma_ref( self, sigma ):
        # CESM2 version verified equality in z coordinates. not necessary here
        if True:
            pass
        else:
            raise Exception('Use this message to inform of any useful exceptions defining sigma_ref.')
        self._sigma_ref = sigma


    def get_zr( self ):
        # Use ref_rho to assign a reference level to all rho(x,y,z)
        # flip rho_ref to make rho the coordinate and z the variable
        flipped_ref = xr.DataArray( \
                        data = self.sigma_ref['z'].values ,
                        coords = {'sigma':self.sigma_ref.values} )

        # Now interpolate to values of sigma. this returns z where sigma_Ref = sigma
        zr = flipped_ref.interp( sigma = self.sigma , kwargs = {'fill_value':'extrapolate' } )
        
        # Set zr = 1 meter wherever it is lower than that (model resolution)
        zr = zr.where( zr < 0 , other = np.nan )
        #nan_mask = ~ np.isnan( zr ); # unable to find zr for some depths
        
        zr = zr.where( zr < -5 , other = - 5 )#.where( nan_mask );
        return zr
    
    def displacement_points( self ):
        # Return array of points between zr and z to help with integrals.

        # determine z-range by stacking z and zr along new dimension
        zr_pos = self.zr.expand_dims( {'zprime' : [0] } )
        zr_to = self.ds['zMid'].expand_dims( {'zprime' : [1] } )
        zrange = np.concatenate( [  zr_pos.values, zr_to.values ] , axis = 0 )
        # reorganize onto an xr object
        zrange = xr.DataArray( data = zrange , 
                              dims = zr_pos.dims, 
                             coords = { 'zprime':[0,1] , 
                                        'nVertLevels' : zr_pos['nVertLevels'] ,
                                        'lat' : zr_pos['lat'] , 
                                        'lon' : zr_pos['lon'] } )
        
        zrange = zrange.interp( zprime = np.arange( 0.05, 0.96, 0.1 )  );
        # zrange is now an array with 10 points inbetween zr and z.
        return zrange

    def boussinesq_PI2( self ):
        # Get Boussinesq expression of APE. 
        # begin with vertical displacement
        DZ = self.ds['zMid'] - self.zr ;
        # z-range over which density differences will be integrated
        zrange = self.displacement_points()

        # Evaluate reference and virtual densities at zrange pressures
        # This creates a new dimension zprime
        sigma_0_eval = self.sigma_ref.interp( z = zrange, 
                                method = 'linear' , 
                                kwargs = {'fill_value':'extrapolate'} );
        sigma_prime = self.sigma - sigma_0_eval

        # now integrate over zprime
        pi2 = 9.81 / 1024 * ( sigma_prime ).mean( 'zprime' ) * DZ
        # eliminate negatives (super small due to rounding errors and limited dz resolution)
        pi2 = pi2.where( pi2 > 0, other = 0 )
        return pi2
    
    def sigma_prime( self ):
        # interpolate sigma ref to zmid of each grid cell
        reference = self.sigma_ref.interp( z = self.ds['zMid'] , 
                                method = 'linear' , 
                                kwargs = { 'fill_value' : 'extrapolate' } ).values
        return self.sigma - reference

    def upwelling_work( self ):
        # approximate w at center of cell
        w = self.ds['vertVelocityTop'] # at top of cell
        dwdz = w.differentiate( 'nVertLevelsP1' ) # diff w/ respect to vertical grid
        w_center = w + 0.5 * dwdz; # 0.5 to put in middle of cell
        # evaluate rho' * g * w
        work = 9.81 * self.sigma_prime() * w_center.values
        return work
    
    def diabatic_tracer( self, variable ):
        total = self.ds[ variable + 'Tend' ] / self.ds['layerThickness']
        advec = self.ds[ variable + 'HorizontalAdvectionTendency'] \
                 + self.ds[ variable + 'VerticalAdvectionTendency']
        return total - advec
    
    def boussinesq_diabatic( self ):
        # Compute diabatic changes to Ea 
        DThetaDt = self.diabatic_tracer( 'temperature' )
        DSDt = self.diabatic_tracer( 'salinity' )
        Theta = self.ds['temperature'].persist()
        Salt = self.ds['salinity'].persist()
        # vertical displacement measures
        dz = self.ds['zMid'] - self.zr
        zrange = self.displacement_points()
        
        
        alpha = gsw.alpha( Salt, Theta, np.abs( zrange ) )
        G_theta = -9.81 * alpha.mean( 'zprime' ) * dz 
        
        beta = gsw.beta( Salt, Theta, np.abs( zrange ) )
        G_theta = -9.81 * beta.mean( 'zprime' ) * dz 
        
        heat_contrib = 1024 * G_theta * DThetaDt
        salt_contrib = 1024 * G_salt * DSDt
        
        return heat_contrib, salt_contrib



class Scotti:
    '''
    Use this class to instantiate analyses following the Scotti and White (2014)
    breakdown of energetic interactions between mean state and turbulent fluctuations.
    '''
    def __init__( self,  ds, mean_ds, sigma_ref ):
        self.avg = mean_ds; 
        self.mean_buoy = buoyancy_from_rho( self.avg['potentialDensity'] )
        self.buoy = buoyancy_from_rho( ds['potentialDensity'] )
        self.ds = ds; 
        # Convert sigma reference to buoyancy reference
        self.buoy_ref = buoyancy_from_rho( sigma_ref );
        self.reference_N2 = - self.buoy_ref.differentiate('z')

    def buoy_ref_on_model_levels( self ):
        # Substitute coordinates on reference buoyancy
        # to show nVertLevels instead of z
        nu_ref = self.buoy_ref.interp( z = self.ds['refZMid'].isel( Time = 0 ) )
        nu_ref = xr.DataArray( data = nu_ref.values, dims = ('nVertLevels') )
        return nu_ref


    def find_zr( self , buoyancy_data ):
        # Use ref_rho to assign a reference level to all rho(x,y,z)
        # flip rho_ref to make rho the coordinate and z the variable
        flipped_ref = xr.DataArray( \
                        data = self.ds['refZMid'].isel( Time = 0 ).values ,
                        coords = {'buoy':self.buoy_ref.values} )

        # Now interpolate to values of sigma. this returns z where sigma_Ref = sigma
        zr = flipped_ref.interp( buoy = buoyancy_data , kwargs = {'fill_value':'extrapolate' } )

        # Set zr = 1 meter wherever it is lower than that (model resolution)
        zr = zr.where( zr < 0 , other = np.nan )
        #nan_mask = ~ np.isnan( zr ); # unable to find zr for some depths

        zr = zr.where( zr < -5 , other = - 5 )#.where( nan_mask );
        return zr

    def buoyancy_displacement_points( self, min_buoy, max_buoy ):
        # Given two limit buoyancies, expand arrays and return
        # an object that transitions between the two buoyancies
        # along a new dimensions. This dimension replaces the
        # integration variable S in Scotti and White (2014).
        DB = max_buoy - min_buoy
        min_buoy = min_buoy.expand_dims( { 'bprime':[0] } )
        max_buoy = max_buoy.expand_dims( { 'bprime':[1] } )
        # Now concatenate onto a single object
        buoyancies = xr.concat( [min_buoy, max_buoy] , dim = 'bprime' )
        bp_values = np.arange( 0.025, 0.98, 0.05 );
        # Integration will require multiplying by max_buoy - min_buoy
        
        return buoyancies.interp( bprime = bp_values ), DB

    def mean_Ea( self ):
        # Calculate Ea stored in mean buoyancy field.
        buoy_points, buoy_range = self.buoyancy_displacement_points(  \
               self.buoy_ref_on_model_levels(), self.mean_buoy )
        # Now interpolate z - z_r( b' ) for b' in buoy_points
        z_r_virtual = self.find_zr(  buoy_points )
        # Integrate along buoy_points
        energy = ( self.avg['zMid'] - z_r_virtual ).mean( 'bprime' ) * buoy_range
        # Get rid of spurious, minuscule negative values 
        return delete_negatives( energy ) 

    def turb_Ea( self ):
        # Calculate Ea stored in turbulent buoyancy field
        buoy_points, buoy_range = self.buoyancy_displacement_points( \
                self.mean_buoy, self.buoy )
        # Interpolate zr( b_mean ) - z_r( b' ) for b' in buoy_points
        zr_of_mean = self.find_zr( self.mean_buoy )
        zr_of_displacements = self.find_zr( buoy_points )
        # Now integrate
        energy = buoy_range * ( zr_of_mean - zr_of_displacements ).mean('bprime')
        # Get rid of spurious, minuscule negative values
        return delete_negatives( energy )
    
    def turbulent_ape_production( self ):
        # Estimate turbulent APE production
        # Compute N2 in reference buoyancy profile, and make buoy_ref the coordinate
        N2_for_mean = xr.DataArray( data = self.reference_N2.values, dims = ('buoy'), 
                                   coords = {'buoy' : self.buoy_ref.values } )
        # Evalute N2 at mean buoyancy values, return array with spatial coordinates
        N2_for_mean = N2_for_mean.interp( buoy = self.mean_buoy ); 

        # Compute spatial gradient of mean buoyancy field
        bgrad = vectors.gradient_vec( self.mean_buoy, self.ds )

        # Prepare to multiply by velocity vector and bprime
        uprime = vectors.uvec( self.ds - self.avg )
        bprime = self.buoy - self.mean_buoy

        # Equalize dimensions to allow dot product
        bgrad = bgrad.rename( {'diff_along' : 'direction' } )
        cross_flux = ( ( bprime * uprime ) * bgrad ).sum( 'direction' )

        return - (1/N2_for_mean) * cross_flux




