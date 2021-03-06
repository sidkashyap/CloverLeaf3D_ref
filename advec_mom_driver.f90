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

!>  @brief Momentum advection driver
!>  @author Wayne Gaudin
!>  @details Invokes the user specified momentum advection kernel.

MODULE advec_mom_driver_module

CONTAINS

SUBROUTINE advec_mom_driver(flopCount,mem,tile,which_vel,direction,sweep_number)

  USE clover_module
  USE advec_mom_kernel_mod

  IMPLICIT NONE

  INTEGER :: tile,which_vel,direction,sweep_number
    INTEGER(KIND=8) :: flopCount,mem


  IF(which_vel.EQ.1)THEN
    CALL advec_mom_kernel(flopCount,                            &
                          mem,                                  &
                          chunk%tiles(tile)%t_xmin,            &
                          chunk%tiles(tile)%t_xmax,              &
                          chunk%tiles(tile)%t_ymin,              &
                          chunk%tiles(tile)%t_ymax,              &
                          chunk%tiles(tile)%t_zmin,              &
                          chunk%tiles(tile)%t_zmax,              &
                          chunk%tiles(tile)%field%xvel1,              &
                          chunk%tiles(tile)%field%mass_flux_x,        &
                          chunk%tiles(tile)%field%vol_flux_x,         &
                          chunk%tiles(tile)%field%mass_flux_y,        &
                          chunk%tiles(tile)%field%vol_flux_y,         &
                          chunk%tiles(tile)%field%mass_flux_z,        &
                          chunk%tiles(tile)%field%vol_flux_z,         &
                          chunk%tiles(tile)%field%volume,             &
                          chunk%tiles(tile)%field%density1,           &
                          chunk%tiles(tile)%field%work_array1,        &
                          chunk%tiles(tile)%field%work_array2,        &
                          chunk%tiles(tile)%field%work_array3,        &
                          chunk%tiles(tile)%field%work_array4,        &
                          chunk%tiles(tile)%field%work_array5,        &
                          chunk%tiles(tile)%field%work_array6,        &
                          chunk%tiles(tile)%field%work_array7,        &
                          chunk%tiles(tile)%field%celldx,             &
                          chunk%tiles(tile)%field%celldy,             &
                          chunk%tiles(tile)%field%celldz,             &
                          advect_x,                               &
                          which_vel,                              &
                          sweep_number,                           &
                          direction,                              &
                          BLOCK%x,                           &
                          BLOCK%y,                           &
                          BLOCK%z)
  ELSEIF(which_vel.EQ.2)THEN
    CALL advec_mom_kernel(flopCount,                            &
                          mem,                                  &
                          chunk%tiles(tile)%t_xmin,            &
                          chunk%tiles(tile)%t_xmax,              &
                          chunk%tiles(tile)%t_ymin,              &
                          chunk%tiles(tile)%t_ymax,              &
                          chunk%tiles(tile)%t_zmin,              &
                          chunk%tiles(tile)%t_zmax,              &
                          chunk%tiles(tile)%field%yvel1,              &
                          chunk%tiles(tile)%field%mass_flux_x,        &
                          chunk%tiles(tile)%field%vol_flux_x,         &
                          chunk%tiles(tile)%field%mass_flux_y,        &
                          chunk%tiles(tile)%field%vol_flux_y,         &
                          chunk%tiles(tile)%field%mass_flux_z,        &
                          chunk%tiles(tile)%field%vol_flux_z,         &
                          chunk%tiles(tile)%field%volume,             &
                          chunk%tiles(tile)%field%density1,           &
                          chunk%tiles(tile)%field%work_array1,        &
                          chunk%tiles(tile)%field%work_array2,        &
                          chunk%tiles(tile)%field%work_array3,        &
                          chunk%tiles(tile)%field%work_array4,        &
                          chunk%tiles(tile)%field%work_array5,        &
                          chunk%tiles(tile)%field%work_array6,        &
                          chunk%tiles(tile)%field%work_array7,        &
                          chunk%tiles(tile)%field%celldx,             &
                          chunk%tiles(tile)%field%celldy,             &
                          chunk%tiles(tile)%field%celldz,             &
                          advect_x,                               &
                          which_vel,                              &
                          sweep_number,                           &
                          direction,                              &
                          BLOCK%x,                           &
                          BLOCK%y,                           &
                          BLOCK%z                          )
  ELSEIF(which_vel.EQ.3)THEN
    CALL advec_mom_kernel(flopCount,                            &
                          mem,                                  &
                          chunk%tiles(tile)%t_xmin,            &
                          chunk%tiles(tile)%t_xmax,              &
                          chunk%tiles(tile)%t_ymin,              &
                          chunk%tiles(tile)%t_ymax,              &
                          chunk%tiles(tile)%t_zmin,              &
                          chunk%tiles(tile)%t_zmax,              &
                          chunk%tiles(tile)%field%zvel1,              &
                          chunk%tiles(tile)%field%mass_flux_x,        &
                          chunk%tiles(tile)%field%vol_flux_x,         &
                          chunk%tiles(tile)%field%mass_flux_y,        &
                          chunk%tiles(tile)%field%vol_flux_y,         &
                          chunk%tiles(tile)%field%mass_flux_z,        &
                          chunk%tiles(tile)%field%vol_flux_z,         &
                          chunk%tiles(tile)%field%volume,             &
                          chunk%tiles(tile)%field%density1,           &
                          chunk%tiles(tile)%field%work_array1,        &
                          chunk%tiles(tile)%field%work_array2,        &
                          chunk%tiles(tile)%field%work_array3,        &
                          chunk%tiles(tile)%field%work_array4,        &
                          chunk%tiles(tile)%field%work_array5,        &
                          chunk%tiles(tile)%field%work_array6,        &
                          chunk%tiles(tile)%field%work_array7,        &
                          chunk%tiles(tile)%field%celldx,             &
                          chunk%tiles(tile)%field%celldy,             &
                          chunk%tiles(tile)%field%celldz,             &
                          advect_x,                               &
                          which_vel,                              &
                          sweep_number,                           &
                          direction,                              &
                          BLOCK%x,                           &
                          BLOCK%y,                           &
                          BLOCK%z                          )
  ENDIF




END SUBROUTINE advec_mom_driver

END MODULE advec_mom_driver_module
