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

  INTEGER :: sweep_number,direction,c

  INTEGER :: xvel,yvel,zvel

  INTEGER :: fields(NUM_FIELDS)
  INTEGER (KINd=8) :: flopCount,tmpFlopCount


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


  IF(profiler_on) kernel_time=timer()
  CALL update_halo(fields,2)
  IF(profiler_on) profiler%halo_exchange=profiler%halo_exchange+(timer()-kernel_time)

  IF(profiler_on) kernel_time=timer()
  DO c=1,chunks_per_task
    tmpFlopCount=0
    CALL advec_cell_driver(tmpFlopCount,c,sweep_number,direction)
    flopCount=flopCount+tmpFlopCount
  ENDDO
  IF(profiler_on) profiler%cell_advection=profiler%cell_advection+(timer()-kernel_time)

  !PRINT *,'BEFORE: cell flop', profiler%advec_cell_flop,flopCount
  IF(profiler_on) profiler%advec_cell_flop=profiler%advec_cell_flop+flopCount

  !PRINT *,'AFTER: cell flop', profiler%advec_cell_flop,flopCount


  flopCount=0
  fields=0
  fields(FIELD_DENSITY1)=1
  fields(FIELD_ENERGY1)=1
  fields(FIELD_XVEL1)=1
  fields(FIELD_YVEL1)=1
  fields(FIELD_ZVEL1)=1
  fields(FIELD_MASS_FLUX_X)=1
  fields(FIELD_MASS_FLUX_Y)=1
  fields(FIELD_MASS_FLUX_Z)=1
  IF(profiler_on) kernel_time=timer()
  CALL update_halo(fields,2)
  IF(profiler_on) profiler%halo_exchange=profiler%halo_exchange+(timer()-kernel_time)

  IF(profiler_on) kernel_time=timer()
  DO c=1,chunks_per_task
    CALL advec_mom_driver(tmpFlopCount,c,xvel,direction,sweep_number) 
    flopCount=flopCount+tmpFlopCount
  ENDDO
  tmpFlopCount=0
  DO c=1,chunks_per_task
    CALL advec_mom_driver(tmpFlopCount,c,yvel,direction,sweep_number) 
    flopCount=flopCount+tmpFlopCount
  ENDDO
  tmpFlopCount=0
  DO c=1,chunks_per_task
    CALL advec_mom_driver(tmpFlopCount,c,zvel,direction,sweep_number) 
    flopCount=flopCount+tmpFlopCount
  ENDDO
  tmpFlopCount=0

  IF(profiler_on) profiler%mom_advection=profiler%mom_advection+(timer()-kernel_time)

  !PRINT *,'BEFORE: mom flop', profiler%advec_mom_flop,flopCount
  IF(profiler_on) profiler%advec_mom_flop=profiler%advec_mom_flop+flopCount

  !PRINT *,'AFTER: mom flop', profiler%advec_mom_flop,flopCount

  flopCount=0
  sweep_number=2
  direction=g_ydir

  fields=0
  fields(FIELD_ENERGY1)=1
  fields(FIELD_DENSITY1)=1
  fields(FIELD_VOL_FLUX_X)=1
  fields(FIELD_VOL_FLUX_Y)=1
  fields(FIELD_VOL_FLUX_Z)=1
  IF(profiler_on) kernel_time=timer()
  CALL update_halo(fields,2)
  IF(profiler_on) profiler%halo_exchange=profiler%halo_exchange+(timer()-kernel_time)
  IF(profiler_on) kernel_time=timer()
  
  flopCount=0
  
  DO c=1,chunks_per_task
    tmpFlopCount=0
    CALL advec_cell_driver(tmpFlopCount,c,sweep_number,direction)
    flopCount=flopCount+tmpFlopCount
  ENDDO
  IF(profiler_on) profiler%cell_advection=profiler%cell_advection+(timer()-kernel_time)

  !PRINT *,'BEFORE: cell flop', profiler%advec_cell_flop,flopCount
  IF(profiler_on) profiler%advec_cell_flop=profiler%advec_cell_flop+flopCount

  !PRINT *,'AFTER: cell flop', profiler%advec_cell_flop,flopCount

  fields=0
  fields(FIELD_DENSITY1)=1
  fields(FIELD_ENERGY1)=1
  fields(FIELD_XVEL1)=1
  fields(FIELD_YVEL1)=1
  fields(FIELD_ZVEL1)=1
  fields(FIELD_MASS_FLUX_X)=1
  fields(FIELD_MASS_FLUX_Y)=1
  fields(FIELD_MASS_FLUX_Z)=1
  IF(profiler_on) kernel_time=timer()
  CALL update_halo(fields,2)
  IF(profiler_on) profiler%halo_exchange=profiler%halo_exchange+(timer()-kernel_time)

  flopCount=0
  IF(profiler_on) kernel_time=timer()
  DO c=1,chunks_per_task
    tmpFlopCount=0
    CALL advec_mom_driver(tmpFlopCount,c,xvel,direction,sweep_number) 
    flopCount=flopCount+tmpFlopCount
  ENDDO


  DO c=1,chunks_per_task

    tmpFlopCount=0
    CALL advec_mom_driver(tmpFlopCount,c,yvel,direction,sweep_number) 
    flopCount=flopCount+tmpFlopCount
  ENDDO


  DO c=1,chunks_per_task
    tmpFlopCount=0 
    CALL advec_mom_driver(tmpFlopCount,c,zvel,direction,sweep_number) 
    flopCount=flopCount+tmpFlopCount
  ENDDO

  IF(profiler_on) profiler%mom_advection=profiler%mom_advection+(timer()-kernel_time)

  !PRINT *,'BEFORE: mom flop', profiler%advec_mom_flop,flopCount
  IF(profiler_on) profiler%advec_mom_flop=profiler%advec_mom_flop+flopCount
  !PRINT *,'AFTER: mom flop', profiler%advec_mom_flop,flopCount

  sweep_number=3
  IF(advect_x)      direction=g_zdir
  IF(.not.advect_x) direction=g_xdir

  flopCount=0
  IF(profiler_on) kernel_time=timer()
  DO c=1,chunks_per_task
    tmpFlopCount=0
    CALL advec_cell_driver(tmpFlopCount,c,sweep_number,direction)
    flopCount=flopCount+tmpFlopCount
  ENDDO
  IF(profiler_on) profiler%cell_advection=profiler%cell_advection+(timer()-kernel_time)

  !PRINT *,'BEFORE: cell flop', profiler%advec_cell_flop,flopCount
  IF(profiler_on) profiler%advec_cell_flop=profiler%advec_cell_flop+flopCount

  !PRINT *,'AFTER: cell flop', profiler%advec_cell_flop,flopCount
  fields=0
  fields(FIELD_DENSITY1)=1
  fields(FIELD_ENERGY1)=1
  fields(FIELD_XVEL1)=1
  fields(FIELD_YVEL1)=1
  fields(FIELD_ZVEL1)=1
  fields(FIELD_MASS_FLUX_X)=1
  fields(FIELD_MASS_FLUX_Y)=1
  fields(FIELD_MASS_FLUX_Z)=1
  IF(profiler_on) kernel_time=timer()
  CALL update_halo(fields,2)
  IF(profiler_on) profiler%halo_exchange=profiler%halo_exchange+(timer()-kernel_time)


  flopCount=0
  tmpFlopCount=0
  IF(profiler_on) kernel_time=timer()
  DO c=1,chunks_per_task
    tmpFlopCount=0
    CALL advec_mom_driver(tmpFlopCount,c,xvel,direction,sweep_number) 
    flopCount=flopCount+tmpFlopCount
  ENDDO

  DO c=1,chunks_per_task
    tmpFlopCount=0
    CALL advec_mom_driver(tmpFlopCount,c,yvel,direction,sweep_number) 
    flopCount=flopCount+tmpFlopCount
  ENDDO
  DO c=1,chunks_per_task
    tmpFlopCount=0
    CALL advec_mom_driver(tmpFlopCount,c,zvel,direction,sweep_number) 
    flopCount=flopCount+tmpFlopCount
  ENDDO

  
  IF(profiler_on) profiler%mom_advection=profiler%mom_advection+(timer()-kernel_time)

  !PRINT *,'BEFORE: mom flop', profiler%advec_mom_flop,flopCount
  IF(profiler_on) profiler%advec_mom_flop=profiler%advec_mom_flop+flopCount

  !PRINT *,'AFTER: mom flop', profiler%advec_mom_flop,flopCount

END SUBROUTINE advection

END MODULE advection_module
