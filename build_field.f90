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

!>  @brief  Allocates the data for each mesh chunk
!>  @author Wayne Gaudin
!>  @details The data fields for the mesh chunk are allocated based on the mesh
!>  size.

MODULE build_field_module

CONTAINS

SUBROUTINE build_field(tile)

   USE clover_module

   IMPLICIT NONE

   INTEGER :: tile
   INTEGER :: x,y,z,bx,by,bz,BLOCK_SIZE_X,BLOCK_SIZE_Y,BLOCK_SIZE_Z
   BLOCK_SIZE_X=BLOCK%x
   BLOCK_SIZE_Y=BLOCK%y
   BLOCK_SIZE_Z=BLOCK%z

   !write(*,*) "Tile ", tile, chunk%tiles(tile)%t_xmin, chunk%tiles(tile)%t_xmax, chunk%tiles(tile)%t_ymin, chunk%tiles(tile)%t_ymax, &
   !         chunk%tiles(tile)%t_zmin, chunk%tiles(tile)%t_zmax



   ALLOCATE(chunk%tiles(tile)%field%density0  (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2, &
                                               chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2, &
                                               chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))
   ALLOCATE(chunk%tiles(tile)%field%density1  (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2, &
                                               chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2, &
                                               chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))
   ALLOCATE(chunk%tiles(tile)%field%energy0   (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2, &
                                               chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2, &
                                               chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))
   ALLOCATE(chunk%tiles(tile)%field%energy1   (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2, &
                                               chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2, &
                                               chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))
   ALLOCATE(chunk%tiles(tile)%field%pressure  (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2, &
                                               chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2, &
                                               chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))
   ALLOCATE(chunk%tiles(tile)%field%viscosity (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2, &
                                               chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2, &
                                               chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))
   ALLOCATE(chunk%tiles(tile)%field%soundspeed(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2, &
                                               chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2, &
                                               chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))

   ALLOCATE(chunk%tiles(tile)%field%xvel0(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                          chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                          chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))
   ALLOCATE(chunk%tiles(tile)%field%xvel1(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                          chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                          chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))
   ALLOCATE(chunk%tiles(tile)%field%yvel0(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                          chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                          chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))
   ALLOCATE(chunk%tiles(tile)%field%yvel1(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                          chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                          chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))
   ALLOCATE(chunk%tiles(tile)%field%zvel0(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                          chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                          chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))
   ALLOCATE(chunk%tiles(tile)%field%zvel1(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                          chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                          chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))


   ALLOCATE(chunk%tiles(tile)%field%vol_flux_x (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                                chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2, &
                                                chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))
   ALLOCATE(chunk%tiles(tile)%field%mass_flux_x(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                                chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2, &
                                                chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))
   ALLOCATE(chunk%tiles(tile)%field%vol_flux_y (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2, &
                                                chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                                chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))
   ALLOCATE(chunk%tiles(tile)%field%mass_flux_y(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2, &
                                                chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                                chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))
   ALLOCATE(chunk%tiles(tile)%field%vol_flux_z (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2, &
                                                chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2, &
                                                chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))
   ALLOCATE(chunk%tiles(tile)%field%mass_flux_z(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2, &
                                                chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2, &
                                                chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))

   ALLOCATE(chunk%tiles(tile)%field%work_array1(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                                chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                                chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))
   ALLOCATE(chunk%tiles(tile)%field%work_array2(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                                chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                                chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))
   ALLOCATE(chunk%tiles(tile)%field%work_array3(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                                chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                                chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))
   ALLOCATE(chunk%tiles(tile)%field%work_array4(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                                chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                                chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))
   ALLOCATE(chunk%tiles(tile)%field%work_array5(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                                chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                                chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))
   ALLOCATE(chunk%tiles(tile)%field%work_array6(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                                chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                                chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))
   ALLOCATE(chunk%tiles(tile)%field%work_array7(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                            chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                            chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))

   ALLOCATE(chunk%tiles(tile)%field%cellx   (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2))
   ALLOCATE(chunk%tiles(tile)%field%celly   (chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2))
   ALLOCATE(chunk%tiles(tile)%field%cellz   (chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))
   ALLOCATE(chunk%tiles(tile)%field%vertexx (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3))
   ALLOCATE(chunk%tiles(tile)%field%vertexy (chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3))
   ALLOCATE(chunk%tiles(tile)%field%vertexz (chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))
   ALLOCATE(chunk%tiles(tile)%field%celldx  (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2))
   ALLOCATE(chunk%tiles(tile)%field%celldy  (chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2))
   ALLOCATE(chunk%tiles(tile)%field%celldz  (chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))
   ALLOCATE(chunk%tiles(tile)%field%vertexdx(chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3))
   ALLOCATE(chunk%tiles(tile)%field%vertexdy(chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3))
   ALLOCATE(chunk%tiles(tile)%field%vertexdz(chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))
   ALLOCATE(chunk%tiles(tile)%field%volume  (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2, &
                                             chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2, &
                                             chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))
   ALLOCATE(chunk%tiles(tile)%field%xarea   (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+3, &
                                             chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2, &
                                             chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))
   ALLOCATE(chunk%tiles(tile)%field%yarea   (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2, &
                                             chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+3, &
                                             chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+2))
   ALLOCATE(chunk%tiles(tile)%field%zarea   (chunk%tiles(tile)%t_xmin-2:chunk%tiles(tile)%t_xmax+2, &
                                             chunk%tiles(tile)%t_ymin-2:chunk%tiles(tile)%t_ymax+2, &
                                             chunk%tiles(tile)%t_zmin-2:chunk%tiles(tile)%t_zmax+3))

   ! Zeroing isn't strictly neccessary but it ensures physical pages
   ! are allocated. This prevents first touch overheads in the main code
   ! cycle which can skew timings in the first step

   !$OMP PARALLEL

!$OMP DO COLLAPSE(3)
DO bz=chunk%tiles(tile)%t_zmin-2, BLOCK_SIZE_Z
DO by=chunk%tiles(tile)%t_ymin-2, BLOCK_SIZE_Y
DO bx=chunk%tiles(tile)%t_xmin-2, BLOCK_SIZE_X

DO z=bz, min(BLOCK_SIZE_Z+bz-1,chunk%tiles(tile)%t_zmax+3)
DO y=by, min(BLOCK_SIZE_y+by-1,chunk%tiles(tile)%t_ymax+3)
DO x=bx, min(BLOCK_SIZE_y+bx-1, chunk%tiles(tile)%t_xmax+3)

     chunk%tiles(tile)%field%work_array1(x,y,z)=0.0
     chunk%tiles(tile)%field%work_array2(x,y,z)=0.0
     chunk%tiles(tile)%field%work_array3(x,y,z)=0.0
     chunk%tiles(tile)%field%work_array4(x,y,z)=0.0
     chunk%tiles(tile)%field%work_array5(x,y,z)=0.0
     chunk%tiles(tile)%field%work_array6(x,y,z)=0.0
     chunk%tiles(tile)%field%work_array7(x,y,z)=0.0
END DO
END DO
END DO
END DO
END DO
END DO
!$OMP END DO

!$OMP DO COLLAPSE(3)

DO bz=chunk%tiles(tile)%t_zmin-2, BLOCK_SIZE_Z
DO by=chunk%tiles(tile)%t_ymin-2, BLOCK_SIZE_Y
DO bx=chunk%tiles(tile)%t_xmin-2, BLOCK_SIZE_X

DO z=bz, min(BLOCK_SIZE_Z+bz-1,chunk%tiles(tile)%t_zmax+2)
DO y=by, min(BLOCK_SIZE_y+by-1,chunk%tiles(tile)%t_ymax+2)
DO x=bx, min(BLOCK_SIZE_y+bx-1, chunk%tiles(tile)%t_xmax+2)


   chunk%tiles(tile)%field%density0(x,y,z)=0.0
   chunk%tiles(tile)%field%density1(x,y,z)=0.0
   chunk%tiles(tile)%field%energy0(x,y,z)=0.0
   chunk%tiles(tile)%field%energy1(x,y,z)=0.0
   chunk%tiles(tile)%field%pressure(x,y,z)=0.0
   chunk%tiles(tile)%field%viscosity(x,y,z)=0.0
   chunk%tiles(tile)%field%soundspeed(x,y,z)=0.0

   chunk%tiles(tile)%field%volume(x,y,z)=0.0
   chunk%tiles(tile)%field%xarea(x,y,z)=0.0
   chunk%tiles(tile)%field%yarea(x,y,z)=0.0


END DO
END DO
END DO
END DO
END DO
END DO

!$OMP END DO

!$OMP DO COLLAPSE(3)

DO bz=chunk%tiles(tile)%t_zmin-2, BLOCK_SIZE_Z
DO by=chunk%tiles(tile)%t_ymin-2, BLOCK_SIZE_Y
DO bx=chunk%tiles(tile)%t_xmin-2, BLOCK_SIZE_X

DO z=bz, min(BLOCK_SIZE_Z+bz-1,chunk%tiles(tile)%t_zmax+3)
DO y=by, min(BLOCK_SIZE_y+by-1,chunk%tiles(tile)%t_ymax+3)
DO x=bx, min(BLOCK_SIZE_y+bx-1, chunk%tiles(tile)%t_xmax+3)


   chunk%tiles(tile)%field%xvel0(x,y,z)=0.0
   chunk%tiles(tile)%field%xvel1(x,y,z)=0.0
   chunk%tiles(tile)%field%yvel0(x,y,z)=0.0
   chunk%tiles(tile)%field%yvel1(x,y,z)=0.0
   chunk%tiles(tile)%field%zvel0(x,y,z)=0.0
   chunk%tiles(tile)%field%zvel1(x,y,z)=0.0

END DO
END DO
END DO
END DO
END DO
END DO


!$OMP END DO


!$OMP DO COLLAPSE(3)



DO bz=chunk%tiles(tile)%t_zmin-2, BLOCK_SIZE_Z
DO by=chunk%tiles(tile)%t_ymin-2, BLOCK_SIZE_Y
DO bx=chunk%tiles(tile)%t_xmin-2, BLOCK_SIZE_X

DO z=bz, min(BLOCK_SIZE_Z+bz-1,chunk%tiles(tile)%t_zmax+2)
DO y=by, min(BLOCK_SIZE_y+by-1,chunk%tiles(tile)%t_ymax+2)
DO x=bx, min(BLOCK_SIZE_y+bx-1, chunk%tiles(tile)%t_xmax+2)


   chunk%tiles(tile)%field%vol_flux_x(x,y,z)=0.0
   chunk%tiles(tile)%field%mass_flux_x(x,y,z)=0.0
   chunk%tiles(tile)%field%vol_flux_y(x,y,z)=0.0
   chunk%tiles(tile)%field%mass_flux_y(x,y,z)=0.0


END DO
END DO
END DO
END DO
END DO
END DO

!$OMP END DO
   
!$OMP DO COLLAPSE(3)

DO bz=chunk%tiles(tile)%t_zmin-2, BLOCK_SIZE_Z
DO by=chunk%tiles(tile)%t_ymin-2, BLOCK_SIZE_Y
DO bx=chunk%tiles(tile)%t_xmin-2, BLOCK_SIZE_X

DO z=bz, min(BLOCK_SIZE_Z+bz-1,chunk%tiles(tile)%t_zmax+3)
DO y=by, min(BLOCK_SIZE_y+by-1,chunk%tiles(tile)%t_ymax+3)
DO x=bx, min(BLOCK_SIZE_y+bx-1, chunk%tiles(tile)%t_xmax+3)


   chunk%tiles(tile)%field%vol_flux_z(x,y,z)=0.0
   chunk%tiles(tile)%field%mass_flux_z(x,y,z)=0.0
   
   chunk%tiles(tile)%field%zarea(x,y,z)=0.0


END DO
END DO
END DO
END DO
END DO
END DO

!$OMP END DO

!$OMP END PARALLEL

   chunk%tiles(tile)%field%cellx=0.0
   chunk%tiles(tile)%field%celly=0.0
   chunk%tiles(tile)%field%cellz=0.0
   chunk%tiles(tile)%field%vertexx=0.0
   chunk%tiles(tile)%field%vertexy=0.0
   chunk%tiles(tile)%field%vertexz=0.0
   chunk%tiles(tile)%field%celldx=0.0
   chunk%tiles(tile)%field%celldy=0.0
   chunk%tiles(tile)%field%celldz=0.0
   chunk%tiles(tile)%field%vertexdx=0.0
   chunk%tiles(tile)%field%vertexdy=0.0
   chunk%tiles(tile)%field%vertexdz=0.0



  
END SUBROUTINE build_field

END MODULE build_field_module
