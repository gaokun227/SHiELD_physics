!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Atmosphere Null Model Component.
!*
!* Atmos Null is free software: you can redistribute it and/or modify it
!* under the terms of the GNU Lesser General Public License as published
!* by the Free Software Foundation, either version 3 of the License, or
!* (at your option) any later version.
!*
!* Atmos Null is distributed in the hope that it will be useful, but
!* WITHOUT ANY WARRANTY; without even the implied warranty of
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!* General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with Atmos Null.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

module atmos_tracer_driver_mod
! <CONTACT EMAIL="William.Cooke@noaa.gov">
!   William Cooke
! </CONTACT>

! <REVIEWER EMAIL="Matthew.Harrison@noaa.gov">
!   Matt Harrison
! </REVIEWER>

! <REVIEWER EMAIL="Bruce.Wyman@noaa.gov">
!   Bruce Wyman
! </REVIEWER>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!     This code allows the user to easily add tracers to the FMS framework.
! </OVERVIEW>

! <DESCRIPTION>
!
!     This is a null module to be used when no atmospheric model is required.
!
!     This code allows a user to easily implement tracer code in the FMS
!     framework.  The tracer and tracer tendency arrays are supplied along
!     with longtitude,  latitude, wind, temperature, and pressure
!     information which allows a  user to implement sources and sinks of the
!     tracer which depend on these parameters.
!
!     In the following example, radon being implemented in the atmosphere
!     will be used as an example of how to implement a tracer in the FMS
!     framework.
!
!     Within the global scope of tracer_driver_mod a use
!     statement should be inserted for each tracer to be added.
!<PRE>      use radon_mod, only : radon_sourcesink, radon_init, radon_end </PRE>
!     
!     An integer parameter, which  will be used as an identifier for the
!     tracer, should be assigned.
!<PRE>
!      integer :: nradon 
!</PRE>
!     Within tracer_driver_init a call to the tracer manager is needed in
!     order  to identify which tracer it has set the tracer as.
!<PRE>
!      nradon = get_tracer_index(MODEL_ATMOS,'radon')
!</PRE>
!     Here MODEL_ATMOS is a parameter defined in field_manager. 
!          'radon' is the name of the tracer within the field_table.
!     
!     If the tracer exists then the integer returned will be positive and it
!     can be used to call the initialization routines for the individual
!     tracers.
!<PRE>
!      if (nradon > 0) then
!           call radon_init(Argument list)
!      endif
!</PRE>     
!
!     Within tracer_driver the user can also use the identifier to surround
!     calls to the source-sink routines for the tracer of interest.
!
!<PRE>
!      if (nradon > 0 .and. nradon <= nt) then
!          call radon_sourcesink (Argument list)
!          rdt(:,:,:,nradon)=rdt(:,:,:,nradon)+rtnd(:,:,:)
!      endif
!</PRE>
!
!     It is the users responsibility to add the tendency generated by the
!     sourcesink routine.
!      
!     Within tracer_driver_end the user can add calls to the
!     terminators for the appropriate source sink routines.
!
!<PRE>      call radon_end</PRE>
!
!     This may simply be a deallocation statement or a routine to send
!     output to the logfile stating that the termination routine has been
!     called.
! </DESCRIPTION>


!-----------------------------------------------------------------------

use platform_mod, only: r8_kind, r4_kind
use     time_manager_mod, only : time_type

implicit none
private
!-----------------------------------------------------------------------
!----- interfaces -------

public  atmos_tracer_driver, atmos_tracer_flux_init
public  atmos_tracer_driver_init, atmos_tracer_driver_end
public  atmos_tracer_driver_gather_data
public  atmos_tracer_driver_gather_data_down
public  atmos_tracer_has_surf_setl_flux, get_atmos_tracer_surf_setl_flux

!-----------------------------------------------------------------------
!----------- namelist -------------------
!-----------------------------------------------------------------------
!
!  When initializing additional tracers, the user needs to make the
!  following changes.
!
!  Add an integer variable below for each additional tracer. 
!  This should be initialized to zero. 
!
!-----------------------------------------------------------------------

character(len=6), parameter :: module_name = 'tracer'

logical :: module_is_initialized = .FALSE.


!---- version number -----
character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'
!-----------------------------------------------------------------------

contains

!#######################################################################

! <SUBROUTINE NAME="atmos_tracer_driver">
!   <OVERVIEW>
!     A routine which allows tracer code to be called.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine calls the source sink routines for atmospheric
!     tracers. This is the interface between the dynamical core of the 
!     model and the tracer code. It should supply all the necessary 
!     information to a user that they need in order to calculate the 
!     tendency of that tracer with respect to emissions or chemical losses.
!
!   </DESCRIPTION>
!   <TEMPLATE>
!     call atmos_tracer_driver (is, ie, js, je, Time, lon, lat, land, phalf, pfull, r,  &
!                           u, v, t, q, u_star, rdt, rm, rdiag, kbot)
!   </TEMPLATE>
!   <IN NAME="is, ie, js, je" TYPE="integer">
!     Local domain boundaries.
!   </IN>
!   <IN NAME="Time" TYPE="type(time_type)">
!     Model time.
!   </IN>
!   <IN NAME="lon" TYPE="real" DIM="(:,:)">
!     Longitude of the centre of the model gridcells
!   </IN>
!   <IN NAME="lat" TYPE="real" DIM="(:,:)">
!     Latitude of the centre of the model gridcells
!   </IN>
!   <IN NAME="land" TYPE="logical" DIM="(:,:)">
!     Land/sea mask.
!   </IN>
!   <IN NAME="phalf" TYPE="real" DIM="(:,:,:)">
!     Pressures on the model half levels.
!   </IN>
!   <IN NAME="pfull" TYPE="real" DIM="(:,:,:)">
!     Pressures on the model full levels.
!   </IN>
!   <IN NAME="r" TYPE="real" DIM="(:,:,:,:)">
!     The tracer array in the component model.
!   </IN>
!   <IN NAME="u" TYPE="real" DIM="(:,:,:)">
!     Zonal wind speed.
!   </IN>
!   <IN NAME="v" TYPE="real" DIM="(:,:,:)">
!     Meridonal wind speed.
!   </IN>
!   <IN NAME="t" TYPE="real" DIM="(:,:,:)">
!     Temperature.
!   </IN>
!   <IN NAME="q" TYPE="real" DIM="(:,:,:)">
!     Specific humidity. This may also be accessible as a
!                        portion of the tracer array.
!   </IN>
!   <IN NAME="u_star" TYPE="real" DIM="(:,:)">
!     Friction velocity :: 
!     The magnitude of the wind stress is density*(ustar**2)
!     The drag coefficient for momentum is u_star**2/(u**2+v**2)
!   </IN>
!   <INOUT NAME="rdt" TYPE="real" DIM="(:,:,:,:)">
!     The tendency of the tracer array in the compenent
!     model. The tendency due to sources and sinks computed
!     in the individual tracer routines should be added to
!     this array before exiting tracer_driver.
!   </INOUT>
!   <IN NAME="rm" TYPE="real" DIM="(:,:,:,:)">
!     The tracer array in the component model for the previous timestep.
!   </IN>
!   <INOUT NAME="rdiag" TYPE="real" DIM="(:,:,:,:)">
!     The array of diagnostic tracers. As these may be changed within the
!     tracer routines for diagnostic purposes, they need to be writable.
!   </INOUT>
!   <IN NAME="kbot" TYPE="integer, optional" DIM="(:,:)">
!     Integer array describing which model layer intercepts the surface.
!   </IN>
 subroutine atmos_tracer_driver (is, ie, js, je, Time, lon, lat,  &
                           area, z_pbl, rough_mom, &
                           land, phalf, pfull,     &
                           u, v, t, q, r,          &
                           rm, rdt, dt,     &
                           u_star, b_star, q_star, &
                           z_half, z_full,&
                           t_surf_rad, albedo, &
                           Time_next, mask, &
                           kbot)

!-----------------------------------------------------------------------
integer, intent(in)                           :: is, ie, js, je
type(time_type), intent(in)                   :: Time
real, intent(in),    dimension(:,:)           :: lon, lat
real, intent(in),    dimension(:,:)           :: u_star, b_star, q_star
real, intent(in),    dimension(:,:)           :: land
real, intent(in),    dimension(:,:)           :: area, z_pbl, rough_mom
real, intent(in),    dimension(:,:,:)         :: phalf, pfull
real, intent(in),    dimension(:,:,:)         :: u, v, t, q
real, intent(inout), dimension(:,:,:,:)       :: r
real, intent(inout), dimension(:,:,:,:)       :: rm
real, intent(inout), dimension(:,:,:,:)       :: rdt
real, intent(in)                              :: dt !timestep(used in chem_interface)
real, intent(in),    dimension(:,:,:)         :: z_half !height in meters at half levels
real, intent(in),    dimension(:,:,:)         :: z_full !height in meters at full levels
real, intent(in),    dimension(:,:)           :: t_surf_rad !surface temperature
real, intent(in),    dimension(:,:)           :: albedo
type(time_type), intent(in)                   :: Time_next
integer, intent(in), dimension(:,:), optional :: kbot
real, intent(in), dimension(:,:,:),  optional :: mask


 return

 end subroutine atmos_tracer_driver
! </SUBROUTINE>

!#######################################################################

! <SUBROUTINE NAME="atmos_tracer_flux_init">
!   <OVERVIEW>
!     Subroutine to initialize the ocean-atmosphere gas flux modules
!   </OVERVIEW>
!   <DESCRIPTION>
!     Subroutine to initialize the ocean-atmosphere gas flux modules
!   </DESCRIPTION>

subroutine atmos_tracer_flux_init

return

end subroutine atmos_tracer_flux_init
! </SUBROUTINE>

!#######################################################################

! <SUBROUTINE NAME="atmos_tracer_driver_init">
!   <OVERVIEW>
!     Subroutine to initialize the tracer driver module.
!   </OVERVIEW>
!   <DESCRIPTION>
!   The purpose of the arguments here are for passing on to the individual
!   tracer code. The user may wish to provide initial values which can be
!   implemented in the initialization part of the tracer code. Remember that
!   the tracer manager will provide a simple fixed or exponential profile if
!   the user provides data for this within the field table. However if a more
!   complicated profile is required then it should be set up in the
!   initialization section of the user tracer code.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call atmos_tracer_driver_init (lonb,latb, r, mask, axes, Time)
!   </TEMPLATE>
!   <IN NAME="lonb" TYPE="real" DIM="(:)">
!     The longitudes for the local domain.
!   </IN>
!   <IN NAME="latb" TYPE="real" DIM="(:)">
!     The latitudes for the local domain.
!   </IN>
!   <IN NAME="mask" TYPE="real, optional" DIM="(:,:,:)">
!      optional mask (0. or 1.) that designates which grid points
!           are above (=1.) or below (=0.) the ground dimensioned as
!           (nlon,nlat,nlev).
!   </IN>
!   <IN NAME="Time" TYPE="type(time_type)">
!     Model time.
!   </IN>
!   <IN NAME="axes" TYPE="integer" DIM="(4)">
!     The axes relating to the tracer array dimensioned as
!      (nlon, nlat, nlev, ntime)
!   </IN>
!   <INOUT NAME="r" TYPE="real" DIM="(:,:,:,:)">
!     Tracer fields dimensioned as (nlon,nlat,nlev,ntrace). 
!   </INOUT>
 subroutine atmos_tracer_driver_init (lonb, latb, r, axes, Time, phalf, mask)

!-----------------------------------------------------------------------
           real, intent(in),    dimension(:,:)             :: lonb, latb
           real, intent(inout), dimension(:,:,:,:)         :: r
type(time_type), intent(in)                                :: Time
        integer, intent(in)                                :: axes(4)
           real, intent(in),    dimension(:,:,:)           :: phalf
           real, intent(in),    dimension(:,:,:), optional :: mask

!-----------------------------------------------------------------------
!
!  When initializing additional tracers, the user needs to make changes 
!
!-----------------------------------------------------------------------

      module_is_initialized = .TRUE.

      return

 end subroutine atmos_tracer_driver_init
! </SUBROUTINE>

!#######################################################################

! <SUBROUTINE NAME="atmos_tracer_driver_end">
!   <OVERVIEW>
!     Subroutine to terminate the tracer driver module.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Termination routine for tracer_driver. It should also call
!     the destructors for the individual tracer routines.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call atmos_tracer_driver_end
!   </TEMPLATE>
 subroutine atmos_tracer_driver_end

!-----------------------------------------------------------------------

      module_is_initialized = .FALSE.

!-----------------------------------------------------------------------

   return

 end subroutine atmos_tracer_driver_end
! </SUBROUTINE>

!#######################################################################

! <SUBROUTINE NAME="atmos_tracer_driver_gather_data">
!   <OVERVIEW>
!     Subroutine to terminate the tracer driver module.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Termination routine for tracer_driver. It should also call
!     the destructors for the individual tracer routines.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call atmos_tracer_driver_gather_data
!   </TEMPLATE>
 subroutine atmos_tracer_driver_gather_data(gas_fields, tr_bot)

use coupler_types_mod, only: coupler_2d_bc_type

type(coupler_2d_bc_type), intent(inout) :: gas_fields
real(r8_kind), dimension(:,:,:), intent(in)      :: tr_bot

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

 return

 end subroutine atmos_tracer_driver_gather_data

!######################################################################
 subroutine atmos_tracer_driver_gather_data_down(gas_fields, tr_bot)

 use coupler_types_mod, only: coupler_2d_bc_type

 type(coupler_2d_bc_type), intent(inout) :: gas_fields
 real(r8_kind), dimension(:,:,:), intent(in)      :: tr_bot

 return

 end subroutine atmos_tracer_driver_gather_data_down
! </SUBROUTINE>
!######################################################################
! given a tracer index, returns true if this tracer has non-zero
! sedimentation flux at the bottom of the atmosphere
function atmos_tracer_has_surf_setl_flux(tr) result(ret)
   logical :: ret
   integer, intent(in) :: tr ! tracer index

   ret=.FALSE.
end function

!######################################################################
! given a tracer index, returns sedimentation flux at the bottom of
! the atmosphere for this tracer
subroutine get_atmos_tracer_surf_setl_flux(tr, setl_flux, dsetl_dtr)
  integer, intent(in)  :: tr         ! tracer index
  real(r8_kind),    intent(out) :: setl_flux(:,:) ! sedimentation flux at the bottom of the atmosphere
  real(r8_kind),    intent(out) :: dsetl_dtr(:,:) ! derivative of sedimentation flux w.r.t.
       ! the tracer concentration in the bottom layer

  setl_flux(:,:) = 0.0 ; dsetl_dtr(:,:) = 0.0
  return
end subroutine
!######################################################################


end module atmos_tracer_driver_mod
