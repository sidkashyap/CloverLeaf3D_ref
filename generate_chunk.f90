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

!>  @brief Mesh chunk generation driver
!>  @author Wayne Gaudin
!>  @details Invoked the users specified chunk generator.

SUBROUTINE generate_chunk(tile)

    USE clover_module
    USE generate_chunk_kernel_module

    IMPLICIT NONE

    INTEGER         :: tile

    INTEGER         :: state
    REAL(KIND=8), DIMENSION(number_of_states) :: state_density,state_energy,state_xvel,state_yvel,state_zvel
    REAL(KIND=8), DIMENSION(number_of_states) :: state_xmin,state_xmax,state_ymin,state_ymax,state_zmin,state_zmax,state_radius
    INTEGER,      DIMENSION(number_of_states) :: state_geometry

    !Sid -> number of states is defined in read_input.f90 -> read the clover.in file to get the state value and choose between the
    !max of that and default state_max
    ! state 1 density=0.2 energy=1.0
    ! state 2 density=1.0 energy=2.5 geometry=cuboid xmin=0.0 xmax=5.0 ymin=0.0 ymax=2.0 zmin=0.0 zmax=4.0
    DO state=1,number_of_states
        state_density(state)=states(state)%density
        state_energy(state)=states(state)%energy
        state_xvel(state)=states(state)%xvel
        state_yvel(state)=states(state)%yvel
        state_zvel(state)=states(state)%zvel
        state_xmin(state)=states(state)%xmin
        state_xmax(state)=states(state)%xmax
        state_ymin(state)=states(state)%ymin
        state_ymax(state)=states(state)%ymax
        state_zmin(state)=states(state)%zmin
        state_zmax(state)=states(state)%zmax
        state_radius(state)=states(state)%radius
        state_geometry(state)=states(state)%geometry
    ENDDO

    CALL generate_chunk_kernel(chunk%tiles(tile)%t_xmin,             &
                               chunk%tiles(tile)%t_xmax,             &
                               chunk%tiles(tile)%t_ymin,             &
                               chunk%tiles(tile)%t_ymax,             &
                               chunk%tiles(tile)%t_zmin,             &
                               chunk%tiles(tile)%t_zmax,             &
                               chunk%tiles(tile)%field%vertexx,           &
                               chunk%tiles(tile)%field%vertexy,           &
                               chunk%tiles(tile)%field%vertexz,           &
                               chunk%tiles(tile)%field%cellx,             &
                               chunk%tiles(tile)%field%celly,             &
                               chunk%tiles(tile)%field%cellz,             &
                               chunk%tiles(tile)%field%density0,          &
                               chunk%tiles(tile)%field%energy0,           &
                               chunk%tiles(tile)%field%xvel0,             &
                               chunk%tiles(tile)%field%yvel0,             &
                               chunk%tiles(tile)%field%zvel0,             &
                               number_of_states,                      &
                               state_density,                         &
                               state_energy,                          &
                               state_xvel,                            &
                               state_yvel,                            &
                               state_zvel,                            &
                               state_xmin,                            &
                               state_xmax,                            &
                               state_ymin,                            &
                               state_ymax,                            &
                               state_zmin,                            &
                               state_zmax,                            &
                               state_radius,                          &
                               state_geometry,                        &
                               g_rect,                                &
                               g_circ,                                &
                               g_point)


END SUBROUTINE generate_chunk
