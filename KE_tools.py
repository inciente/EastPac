import sys, gsw 
import matplotlib.pyplot as plt; 
import numpy as np; 
import pandas as pd; import xarray as xr; 
sys.path.append('/home/noelgb/repositories/EastPac/')
import intercomparison, thermo
sys.path.append('/home/noelgb/notebooks/config_files')
from open_jsons import *
from euc_tools import *
# --------------

''' 
Functions to compute the energy budget of the upper Pacific Ocean
'''

def meridional_KE_flux( pop_ds, integrate = True ):
    # computes KE advective flux in given direction
    # pop_ds is xr.dataset with pop2 output
    # direction determines advective variable
    direction = 'VVEL'
    KE = 1025 / 2 * ( pop_ds['UVEL'] ** 2 + pop_ds['VVEL'] ** 2 \
           + pop_ds['WVEL'].values ** 2 ) * pop_ds[ direction ] / 1e6
    # division over 1e6 adjusts to m/s
    if integrate:
        KE = xz_integral( KE, pop_ds )
    return KE

def meridional_pwork( pop_ds, ref_rho, ref_ssh, integrate = True  ):
    # Compute pressure the meridional pressure work
    vel = pop_ds[ 'VVEL' ] / 100
    # specify reference density and ssh 
    p_anom = thermo.pressure( pop_ds , ref_rho = ref_rho , 
               ssh_val = ref_ssh )
    pwork = vel * p_anom
    if integrate:
        pwork = xz_integral( pwork, pop_ds )
    return pwork

def meridional_APE_flux( pop_ds, ref_rho, integrate = True ):
    # compute APE and its meridional flux
    rho_prime = pop_ds['RHO'] * 1000 - ref_rho; 
    N2_ref = 9.81 / 1025 * ref_rho.differentiate( 'z_t' )
    APE = 9.81 ** 2 * rho_prime ** 2 
    APE = APE / ( 2 * 1025 * N2_ref )
    APE_flux = pop_ds['VVEL'] / 100 * APE
    
    if integrate:
        APE_flux = xz_integral( APE_flux, pop_ds )
    return APE_flux


def wind_work( pop_ds, integrate = True ):
    # wind work on the ocean surface
    power = pop_ds['TAUX'] * pop_ds['UVEL'].isel( z_t = 0 )
    power += pop_ds['TAUY'] * pop_ds['VVEL'].isel( z_t = 0 )
    # now convert to W m-2
    power = power * 1e-7 * 1e4
    if integrate:
        power = xy_integral( power, pop_ds )
    return power  
 

def net_KE( ds, integrate = True ):
    # compute net KE using UVEL2 variables
    ke = ds['UVEL'] + ds['VVEL2'] + ds['WVEL2'].values
    ke = ke / 1e4 * 1025 / 2 
    if integrate:
        ke = xyz_integral( ke, ds )
    return ke

# -------- functions for spatial subsetting and section management

def basin_mask():
    # Load the WOA basin mask at 100 m depth
    path = '/home/noelgb/data/ocean_mask/WOA_basins_100m.nc'
    mask = xr.open_dataset( path ); 
    # change longitude from -180 to 0-360
    nulon = mask['lon'].values;
    negative = nulon < 0; 
    nulon[negative] = nulon[negative] + 360;
    mask['lon'] = nulon;
    # sort 
    mask = mask.sortby('lon')
    return mask

def pacific_only( xr_obj ):
    # Apply WOA mask to object, leaving only data inside the Pacific
    mask = basin_mask(); 
    mask = mask['basins'].interp( lat = xr_obj['lat'] ).interp( \
                        lon = xr_obj['lon'] )
    return xr_obj.where( mask == 2 )

# -------------- CALCULUS TOOLS

def xy_integral( xr_obj, parent_ds ):
    # return area integral (lon, lat)
    data = ( xr_obj * parent_ds['dx'] * parent_ds['dy'] )
    data = data.sum(['lon','lat'])
    return data

def xz_integral( xr_obj, parent_ds ):
    # integrate over lon, depth
    data = ( xr_obj * parent_ds['dx'] * parent_ds['dz'] )
    data = data.sum( ['lon','z_t'] )
    return data

def z_integral( xr_obj, parent_ds ):
    data = ( xr_obj * parent_ds['dz'] ).sum('z_t')
    return data

def xyz_integral( xr_obj , parent_ds ):
    data = xy_integral( xr_obj, parent_ds )
    data = z_integral( data, parent_ds )
    return data 

def horizontal_gradient( xr_obj ):
    lon_deriv = xr_obj.differentiate('lon') / ( \
         110e3 * np.cos( np.pi * xr_obj['lat'] / 180 ) )
    lat_deriv = xr_obj.differentiate('lat') / 110e3
    return lon_deriv, lat_deriv



''' 
NOW ACTUALLY LAY OUT THE COMPUTATION WE'RE INTERESTED IN
'''

def prepare_model_for_energy( model ):
    ds = pull_combined_json( model )
    latitudes = [-35,-15,-3,3,15,35];
    # get metrics and spatial subsetting
    ds = ds.sel( lat = slice( -40, 40 ) );
    ds = ds.sel( lon = slice( 100, 295 ) );
    ds = pacific_only( ds );
    ds = ds.isel( time = range(60) ); # for speed up in testing
    # get reference density and ssh
    ds['rho_ref'] = thermo.reference_density( ds )
    ds['rho_ref'] = ds['rho_ref'].interp( z_t = ds['z_t'] , method = 'linear', 
                                         kwargs = {'fill_value': 'extrapolate' } )
    ds['mean_ssh'] = ds['SSH'].mean(['lon','lat','time']).persist()
    # focus only on the upper 500 m 
    #ds = vertical_chunk( ds , 500 )
    return ds 

def all_merid_fluxes( pop_ds ):
    # compute all meridional fluxes from functions in this module
    fluxes = xr.Dataset()
    fluxes['adv'] = meridional_KE_flux( pop_ds , integrate = True );
    fluxes['pwork'] = meridional_pwork( pop_ds , 
                         ref_rho = pop_ds['rho_ref'], 
                         ref_ssh = pop_ds['mean_ssh'],
                         integrate = True )
    fluxes['APE_adv'] = meridional_APE_flux( pop_ds , 
                         ref_rho = pop_ds['rho_ref'], 
                         integrate = True )
    return fluxes 

def all_conversions( pop_ds ):
    # compute all energy conversions involving area/volume integrals
    # each function returns a time series
    conversions = xr.Dataset()
    conversions['upwelling'] = upwelling_work( pop_ds , 
                                    ref_rho = pop_ds['rho_ref'],
                                    integrate = True )
    conversions['windwork'] = wind_work( pop_ds , integrate = True )
    conversions['mixing'] = mixing_work( pop_ds, integrate = True )
    return conversions

def eval_KE( ds ):
    #ds = prepare_model_for_energy( model )
    latitudes = [-35, -15, -3, 3, 15, 35 ]
    # Get all the components of the energy budget
     
    # ---------- begin with meridional fluxes
    ds_merid = ds.sel( lat = latitudes , method = 'nearest' );
    fluxes = all_merid_fluxes( ds_merid ); 
    
    # ---------- now slice sections where area or volume ints involved
    area_quants = []
    for jj in range( len( latitudes ) - 1 ):
        lat_range = slice( latitudes[jj] , latitudes[jj+1] )
        ds_here = ds.sel( lat = lat_range ); # this will be integrated
        conversions = all_conversions( ds_here )
        conversions = conversions.assign_coords( { 'lat': [latitudes[jj]] } )
        area_quants.append( conversions )
    # concatenate all area and volume integrals
    conversions = xr.concat( area_quants, dim = 'lat' )
    return fluxes, conversions 


''' 
by now (April 15, 2024) it looks like i've settled on working under the Tailleux (2018) formalism.
Therefore, functions above might be deleted at some point. Whatever the case is, the class below 
is becoming the main focus of development now.
'''    


def equality_check( arr1, arr2 ):
    # check that two arrays are equal, element by element

    # first compare lengths
    if ( len( arr1 ) != len( arr2) ):
        return False

    # now check elements
    for i in range( len( arr1 ) ):
      
        if ( arr1[i] != arr2[i] ):
            # Fail one comparison fail whole thing
            return False

    # if you made it this far it means arrays are equal
    return True


class Tailleux:

    '''
Methods designed to compute most variables needed under the Tailleux formalism. 

Input:
    ds    -    snapshot of POP2 dataset (no time dependence) 
    rho_ref    reference density profile (may be computed from thermo.reference_density() )

All quantities will be computed for all spatial domain in ds, so subset ds before instantiating this class.
    '''
    def __init__( self, ds, rho_ref ):
        
        self.ds = ds; 
        self.rho_ref = rho_ref;
        # save rho values, since they will be used again and again
        self.rho = ( self.ds['RHO'] * 1000 ).persist(); 
        self.nan_mask = ~ np.isnan( self.rho ); # useful to have
        # compute zr, since it will also be needed repeatedly
        self.zr = self.get_zr().persist()

    @property
    def ds(self):
        return self._ds

    @ds.setter
    def ds(self, data):
        if 'time' in list( data.dims ):
            raise Exception('Cannot instantiate class. Dataset has time dependence.')
        self._ds = data

    @property
    def rho_ref( self ):
        return self._rho_ref

    @rho_ref.setter
    def rho_ref( self, rho ):
        # make sure that vertical coordinate is same as in ds
        if equality_check( rho['z_t'] , self.ds['z_t'] ):
            pass
        else: 
            raise Exception('Dataset and reference density have different vertical coordinates')
        self._rho_ref = rho


    def get_zr( self ):
        # Use ref_rho to assign a reference level to all rho(x,y,z)
        # flip rho_ref to make rho the coordinate and z the variable
        flipped_ref = xr.DataArray( \
                        data = self.rho_ref['z_t'].values , 
                        coords = {'rho':self.rho_ref.values} )

        # Now interpolate to values of rho
        zr = flipped_ref.interp( rho = self.rho, kwargs = {'fill_value':'extrapolate' } )
        # Set zr = 5 wherever it is lower than that
        nan_mask = ~ np.isnan( zr );
        zr = zr.where( zr > 5 , other = 5 ).where( nan_mask );
        return zr

    def reference_pressure( self ):
        # Compute pressure in reference density profile
        
        # Create xr dataset with dz and rho (so it can be passed to thermo.pressure)
        dummy_data = xr.Dataset() 
        dummy_data['RHO'] = self.rho_ref / 1000; 
        dummy_data['dz'] = thermo.get_dz( self.ds )
        dummy_data['SSH'] = self.ds['SSH'].mean(['lon','lat'])
        
        # compute reference pressure in Pa
        ref_p = thermo.pressure( dummy_data )
        
        return ref_p


    def enthalpy_diff( self ):
        # Compute enthalpy at reference pressure at current AND reference levels
        ref_p = self.reference_pressure() / 1e4 ; # scale to dbar
        p_zr = ref_p.interp( z_t = self.zr ); # array with zr values throughout ds coords 
 
        # Use gsw to get enthalpy with SALT, TEMP, PRESSURE
        ent_simp = gsw.enthalpy_CT_exact( self.ds['SALT'] , self.ds['TEMP'] , ref_p )
        ent_zr =  gsw.enthalpy_CT_exact( self.ds['SALT'] , self.ds['TEMP'] , p_zr );

        return ent_simp - ent_zr 

    def get_PI2( self ):
        # Get full expression of PI2
        # order of zr and z_t reversed because z is defined as depth here
        pi2 = 9.81 * ( self.zr - self.ds['z_t'] ) + self.enthalpy_diff()
        return pi2

    def boussinesq_PI2( self ):
        # Get Boussinesq expression of APE. 
        # begin with vertical displacement
        DZ = self.zr - self.ds['z_t'] ;

        # z-range over which density differences will be integrated
        zrange = xr.concat( [ self.zr.expand_dims( { 'zprime' : [0] } ) , 
                        self.ds['z_t'].expand_dims( {'zprime' : [1] } ) ] , 
                        dim = 'zprime' )
        zrange = zrange.interp( zprime = np.arange( 0.05, 0.96, 0.1 ) ) ;

        # Evaluate reference and virtual densities at zrange pressures
        # This creates a new dimension zprime
        rho_0_eval = self.rho_ref.interp( z_t = zrange ); 
        rho_actual_eval = gsw.rho( self.ds['SALT'] , self.ds['TEMP'] , zrange )
 
        # now integrate over zprime
        pi2 = 9.81 * ( 1 - rho_0_eval / rho_actual_eval ).mean( 'zprime' ) * DZ 
          
        return pi2
        


    def QG_APE( self ):
        # Compute APE using the QG approximation, for comparison purposes
        ape = 9.81 / 2 * ( self.rho - self.rho_ref ) * ( self.zr - self.rho['z_t'] )
        return ape 

    def rho_prime( self ):

        # modified density rho * ( 1 - rho_0 / rho_h ) used in framework
        rho_h = self.rho_ref.interp( z_t = self.zr )

        return self.rho * ( 1 - self.rho_ref / rho_h )

    def upwelling_work( self ):
        # evaluate g * rho' * w 
        work = 9.81 * self.rho_prime() * self.ds['WVEL'].values / 100
        return work   
        
    def advective_term( self, xr_obj, three_d = True ):
        # compute gradient of xr_obj and dot product with advection
        ox, oy = horizontal_gradient( xr_obj )
        advec = self.ds['UVEL'] * ox + self.ds['VVEL'] * oy
        
        if three_d:
            oz = - xr_obj.differentiate('z_t')
            advec = advec + oz * self.ds['WVEL'].values
        
        return advec/100 

    def ML_diabatic( self ):
        # Compute the diabatic change to APE in the mixed layer
        # Get ML mask plus 5 meters
        h = self.ds['HMXL'] / 100 + 5
        ml_mask = self.ds['z_t'] < h
        
        # actual and reference depths
        actual_z = self.ds['z_t'].where( ml_mask ).persist()
        virtual_z = self.zr.where( ml_mask ).persist()

        # Get Theta and Salinity in mixed layer
        theta = self.ds['TEMP'].where( ml_mask ).persist()
        salt = self.ds['SALT'].where( ml_mask ).persist()

        # Compute in-situ temperature and chem potential at z and z_r
        temp_insitu = gsw.t_from_CT( salt, theta, actual_z )
        temp_virtual = gsw.t_from_CT( salt, theta, virtual_z )

        mu_insitu = gsw.chem_potential_water_t_exact( salt, theta, actual_z )
        mu_virtual = gsw.chem_potential_water_t_exact( salt, theta, virtual_z )

        # Compute heat input and salinity changes to mixed layer
        ent_change = self.ds['SHF'] / ( 273.15 + temp_insitu ) / h
        salt_change = - salt * ( np.abs( self.ds['EVAP_F'] ) \
                       - self.ds['PREC_F'] ) / h

        heat_contribution = ( temp_insitu - temp_virtual ) * ent_change
        salt_contribution = ( mu_insitu - mu_virtual ) * salt_change
        
        return heat_contribution, salt_contribution


        # first derivatives of specific entropy
        #det_ds, det_dt = gsw.entropy_first_derivatives( ml_salt.compute() , ml_temp.compute() )












