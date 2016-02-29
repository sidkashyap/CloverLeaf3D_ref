!Crown Copyright 2012 AWE.
!
! This file is part of CloverLeaf.
!
! CloverLeaf is free software: you can redistribute it and/or modify it under 
! the terms of the GNU General Public License as published by the 
! Free Software Foundation, either version 3 of the License, or (at your option) 
! any later version.
!
! CloverLeaf is distributed in the hope that it will be useful, but 
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
! details.
!
! You should have received a copy of the GNU General Public License along with 
! CloverLeaf. If not, see http://www.gnu.org/licenses/.

!>  @brief Top level advection driver
!>  @author Wayne Gaudin
!>  @details Controls the advection step and invokes required communications.

MODULE advection_module

CONTAINS

SUBROUTINE advection()
  !Sid - this is the kernel that is called from hydro.f90

  USE clover_module
  USE advec_cell_driver_module
  USE advec_mom_driver_module
  USE update_halo_module

  IMPLICIT NONE

  INTEGER :: sweep_number,direction,tile

  INTEGER :: xvel,yvel,zvel

  INTEGER :: fields(NUM_FIELDS)
  INTEGER (KIND=8) :: flopCount,tmpFlopCount,mem,tmpMem,tmpFlopCountx,tmpFlopCounty,tmpFlopCountz,tmpMemx,tmpMemy,tmpMemz


  REAL(KIND=8) :: kernel_time,timer

  !Sid - TODO - what is sweep_number?
  sweep_number=1

  !Sid - okay, the direction is either x or z

  IF(advect_x)      direction=g_xdir
  IF(.not.advect_x) direction=g_zdir

  !Sid - initialize the x,y and z velocities
  xvel=g_xdir
  yvel=g_ydir
  zvel=g_zdir

  !Sid - initialize fields 
  fields=0
  fields(FIELD_ENERGY1)=1
  fields(FIELD_DENSITY1)=1
  fields(FIELD_VOL_FLUX_X)=1
  fields(FIELD_VOL_FLUX_Y)=1
  fields(FIELD_VOL_FLUX_Z)=1

  !Sid - variables used to count flop, returned from the kernel
  flopCount=0
  tmpFlopCount=0  
  tmpFlopCountx=0
  tmpFlopCounty=0
  tmpFlopCountz=0
  tmpMemx=0
  tmpMemy=0
  tmpMemz=0
  CALL update_halo(fields,2)


  IF(profiler_on) kernel_time=timer()

  DO tile=1,tiles_per_chunk
    CALL advec_cell_driver(tile,sweep_number,direction)
  ENDDO

  IF(profiler_on) profiler%cell_advection=profiler%cell_advection+(timer()-kernel_time)

  fields=0
  fields(FIELD_DENSITY1)=1
  fields(FIELD_ENERGY1)=1
  fields(FIELD_XVEL1)=1
  fields(FIELD_YVEL1)=1
  fields(FIELD_ZVEL1)=1
  fields(FIELD_MASS_FLUX_X)=1
  fields(FIELD_MASS_FLUX_Y)=1
  fields(FIELD_MASS_FLUX_Z)=1

  CALL update_halo(fields,2)

  IF(profiler_on) kernel_time=timer()
  DO tile=1,tiles_per_chunk
    CALL advec_mom_driver(tmpFlopCountx,tmpMemx,tile,xvel,direction,sweep_number)
    CALL advec_mom_driver(tmpFlopCounty,tmpMemy,tile,yvel,direction,sweep_number)
    CALL advec_mom_driver(tmpFlopCountz,tmpMemz,tile,zvel,direction,sweep_number) 
    flopCount=flopCount+tmpFlopCountx+tmpFlopCounty+tmpFlopCountz
    mem=mem+tmpMemx+tmpMemy+tmpMemz
    tmpFlopCountx=0
    tmpFlopCounty=0
    tmpFlopCountz=0
    tmpMemx=0
    tmpMemy=0
    tmpMemz=0
  ENDDO

  IF(profiler_on) profiler%mom_advection=profiler%mom_advection+(timer()-kernel_time)
  IF(profiler_on) profiler%advec_mom_flop=profiler%advec_mom_flop+flopCount
  IF(profiler_on) profiler%advec_mem=profiler%advec_mem+mem

  !print *,flopCount,":FC ",mem,":mem ",profiler%advec_mom_flop,":tFC ",profiler%advec_mem,":tMem"
  flopCount=0
  mem=0

  sweep_number=2
  direction=g_ydir

  fields=0
  fields(FIELD_ENERGY1)=1
  fields(FIELD_DENSITY1)=1
  fields(FIELD_VOL_FLUX_X)=1
  fields(FIELD_VOL_FLUX_Y)=1
  fields(FIELD_VOL_FLUX_Z)=1

  CALL update_halo(fields,2)


  IF(profiler_on) kernel_time=timer()

  DO tile=1,tiles_per_chunk
    CALL advec_cell_driver(tile,sweep_number,direction)
  ENDDO

  IF(profiler_on) profiler%cell_advection=profiler%cell_advection+(timer()-kernel_time)

  fields=0
  fields(FIELD_DENSITY1)=1
  fields(FIELD_ENERGY1)=1
  fields(FIELD_XVEL1)=1
  fields(FIELD_YVEL1)=1
  fields(FIELD_ZVEL1)=1
  fields(FIELD_MASS_FLUX_X)=1
  fields(FIELD_MASS_FLUX_Y)=1
  fields(FIELD_MASS_FLUX_Z)=1

  CALL update_halo(fields,2)

  IF(profiler_on) kernel_time=timer()
  DO tile=1,tiles_per_chunk


    CALL advec_mom_driver(tmpFlopCountx,tmpMemx,tile,xvel,direction,sweep_number)
    CALL advec_mom_driver(tmpFlopCounty,tmpMemy,tile,yvel,direction,sweep_number)
    CALL advec_mom_driver(tmpFlopCountz,tmpMemz,tile,zvel,direction,sweep_number) 
    flopCount=flopCount+tmpFlopCountx+tmpFlopCounty+tmpFlopCountz
    mem=mem+tmpMemx+tmpMemy+tmpMemz
    tmpFlopCountx=0
    tmpFlopCounty=0
    tmpFlopCountz=0
    tmpMemx=0
    tmpMemy=0
    tmpMemz=0
  ENDDO


  IF(profiler_on) profiler%mom_advection=profiler%mom_advection+(timer()-kernel_time)
  IF(profiler_on) profiler%advec_mom_flop=profiler%advec_mom_flop+flopCount
  IF(profiler_on) profiler%advec_mem=profiler%advec_mem+mem


  !print *,flopCount,":FC ",mem,":mem ",profiler%advec_mom_flop,":tFC ",profiler%advec_mem,":tMem"
  flopCount=0
  mem=0

  !PRINT *,'AFTER: mom flop', profiler%advec_mom_flop,flopCount

  sweep_number=3
  IF(advect_x)      direction=g_zdir
  IF(.not.advect_x) direction=g_xdir

  flopCount=0
  IF(profiler_on) kernel_time=timer()

  DO tile=1,tiles_per_chunk
    CALL advec_cell_driver(tile,sweep_number,direction)
  ENDDO

  IF(profiler_on) profiler%cell_advection=profiler%cell_advection+(timer()-kernel_time)

  fields=0
  fields(FIELD_DENSITY1)=1
  fields(FIELD_ENERGY1)=1
  fields(FIELD_XVEL1)=1
  fields(FIELD_YVEL1)=1
  fields(FIELD_ZVEL1)=1
  fields(FIELD_MASS_FLUX_X)=1
  fields(FIELD_MASS_FLUX_Y)=1
  fields(FIELD_MASS_FLUX_Z)=1

  CALL update_halo(fields,2)


  IF(profiler_on) kernel_time=timer()



  DO tile=1,tiles_per_chunk
    CALL advec_mom_driver(tmpFlopCountx,tmpMemx,tile,xvel,direction,sweep_number)
    CALL advec_mom_driver(tmpFlopCounty,tmpMemy,tile,yvel,direction,sweep_number)
    CALL advec_mom_driver(tmpFlopCountz,tmpMemz,tile,zvel,direction,sweep_number) 
    flopCount=flopCount+tmpFlopCountx+tmpFlopCounty+tmpFlopCountz
    mem=mem+tmpMemx+tmpMemy+tmpMemz
    tmpFlopCountx=0
    tmpFlopCounty=0
    tmpFlopCountz=0
    tmpMemx=0
    tmpMemy=0
    tmpMemz=0
  ENDDO


  IF(profiler_on) profiler%mom_advection=profiler%mom_advection+(timer()-kernel_time)
  IF(profiler_on) profiler%advec_mom_flop=profiler%advec_mom_flop+flopCount
  IF(profiler_on) profiler%advec_mem=profiler%advec_mem+mem
  !print *,flopCount,":FC ",mem,":mem ",profiler%advec_mom_flop,":tFC ",profiler%advec_mem,":tMem"
END SUBROUTINE advection

END MODULE advection_module
