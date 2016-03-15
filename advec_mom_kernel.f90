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

!>  @brief Fortran momentum advection kernel
!>  @author Wayne Gaudin
!>  @details Performs a second order advective remap on the vertex momentum
!>  using van-Leer limiting and directional splitting.
!>  Note that although pre_vol is only set and not used in the update, please
!>  leave it in the method.

MODULE advec_mom_kernel_mod
USE OMP_LIB
CONTAINS

SUBROUTINE advec_mom_kernel(flop,                           &
                            mem,                            &
                            x_min,x_max,y_min,y_max,z_min,z_max, &
                            vel1,                                &
                            mass_flux_x,                         &
                            vol_flux_x,                          &
                            mass_flux_y,                         &
                            vol_flux_y,                          &
                            mass_flux_z,                         &
                            vol_flux_z,                          &
                            volume,                              &
                            density1,                            &
                            node_flux,                           &
                            node_mass_post,                      &
                            node_mass_pre,                       &
                            advec_vel,                           &
                            mom_flux,                            &
                            pre_vol,                             &
                            post_vol,                            &
                            celldx,                              &
                            celldy,                              &
                            celldz,                              &
                            advect_x,                            &
                            which_vel,                           &
                            sweep_number,                        &
                            direction, &
                            BLOCK_SIZE_x,                        &
                            BLOCK_SIZE_y,                        &
                            BLOCK_SIZE_z)

  IMPLICIT NONE
  
  INTEGER :: x_min,x_max,y_min,y_max,z_min,z_max
  INTEGER :: which_vel,sweep_number,direction
  LOGICAL :: advect_x
  INTEGER (KIND=8) :: flop,tmpFlop,mem,tmpMem,iter,limiter_count
  REAL(KIND=8), TARGET,DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3) :: vel1
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2,z_min-2:z_max+2) :: mass_flux_x
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2,z_min-2:z_max+2) :: vol_flux_x
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3,z_min-2:z_max+2) :: mass_flux_y
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3,z_min-2:z_max+2) :: vol_flux_y
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+3) :: mass_flux_z
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+3) :: vol_flux_z
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: volume
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2,z_min-2:z_max+2) :: density1
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3) :: node_flux
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3) :: node_mass_post
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3) :: node_mass_pre
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3) :: advec_vel
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3) :: mom_flux
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3) :: pre_vol
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3,z_min-2:z_max+3) :: post_vol

  REAL(KIND=8), DIMENSION(x_min-2:x_max+2) :: celldx
  REAL(KIND=8), DIMENSION(y_min-2:y_max+2) :: celldy
  REAL(KIND=8), DIMENSION(z_min-2:z_max+2) :: celldz
 
  INTEGER :: j,k,l,bj,bk,bl,nb_x,nb_y,nb_z,block_max,BLOCK_SIZE_x,BLOCK_SIZE_y,BLOCK_SIZE_z
  INTEGER :: upwind,donor,downwind,dif
  REAL(KIND=8) :: sigma,wind,width,tmp1,tmp2,tmp3,tmp4,tmp5
  REAL(KIND=8) :: vdiffuw,vdiffdw,auw,adw,limiter
  REAL(KIND=8) :: vdiffuw2,vdiffdw2,auw2,limiter2, advec_vel_s

  ! Choose the correct velocity, ideally, remove this pointer
  !  if it affects performance.
  ! Leave this one in as a test of performance

  flop=0
  tmpFlop=0

! I think these only have to be done once per cell advection sweep. So put in some logic so they are just done the first time


nb_x=x_max+2
nb_y=y_max+2
nb_z=z_max+2

!$OMP PARALLEL
  IF(sweep_number.EQ.1.AND.direction.EQ.1)THEN ! x first
!$OMP DO collapse(3)

DO bl=z_min-2,nb_z,BLOCK_SIZE_z
    DO bk=y_min-2,nb_y,BLOCK_SIZE_y
      DO bj=x_min-2,nb_x,BLOCK_SIZE_x
        DO l=bl,min(nb_z,BLOCK_SIZE_z+bl-1)
            DO k=bk,min(nb_y,BLOCK_SIZE_y+bk-1)
                DO j=bj,min(nb_x,BLOCK_SIZE_x+bj-1)
                      post_vol(j,k,l)= volume(j,k,l)+vol_flux_y(j  ,k+1,l  )-vol_flux_y(j,k,l) &
                                        +vol_flux_z(j  ,k  ,l+1)-vol_flux_z(j,k,l)
                      pre_vol(j,k,l)=post_vol(j,k,l)+vol_flux_x(j+1,k  ,l  )-vol_flux_x(j,k,l)
                      !print *,"pre_vol(,",j,",",k,",",l,"):",pre_vol(j,k,l)
                    ENDDO
                ENDDO
            ENDDO
        ENDDO
    ENDDO
ENDDO

!$OMP END DO


  iter=(z_max-z_min+5)*(y_max-y_min+5)*(x_max-x_min+5)
  tmpFlop=iter*6 !sc - 16464
  tmpMem=iter*10



  ELSEIF(sweep_number.EQ.1.AND.direction.EQ.3)THEN ! z first
!$OMP DO collapse(3)
DO bl=z_min-2,nb_z,BLOCK_SIZE_z
    DO bk=y_min-2,nb_y,BLOCK_SIZE_y
      DO bj=x_min-2,nb_x,BLOCK_SIZE_x
        DO l=bl,min(nb_z,BLOCK_SIZE_z+bl-1)
            DO k=bk,min(nb_y,BLOCK_SIZE_y+bk-1)
                DO j=bj,min(nb_x,BLOCK_SIZE_x+bj-1)
                      post_vol(j,k,l)= volume(j,k,l)+vol_flux_x(j+1,k  ,l  )-vol_flux_x(j,k,l) &
                                        +vol_flux_y(j  ,k+1,l  )-vol_flux_y(j,k,l)
                      pre_vol(j,k,l)=post_vol(j,k,l)+vol_flux_z(j  ,k  ,l+1)-vol_flux_z(j,k,l)
                        ENDDO

                ENDDO
            ENDDO
        ENDDO
    ENDDO
ENDDO



!$OMP END DO


  iter=(z_max-z_min+1+4)*(y_max-y_min+1+4)*(x_max-x_min+1+4)
  tmpFlop=iter*6
  tmpMem=iter*10


ELSEIF(sweep_number.EQ.2.AND.advect_x)THEN ! x first
!$OMP DO collapse(3)
DO bl=z_min-2,nb_z,BLOCK_SIZE_z
    DO bk=y_min-2,nb_y,BLOCK_SIZE_y

      DO bj=x_min-2,nb_x,BLOCK_SIZE_x
        DO l=bl,min(nb_z,BLOCK_SIZE_z+bl-1)
            DO k=bk,min(nb_y,BLOCK_SIZE_y+bk-1)
                DO j=bj,min(nb_x,BLOCK_SIZE_x+bj-1)
                        post_vol(j,k,l)=volume(j,k,l) +vol_flux_z(j  ,k  ,l+1)-vol_flux_z(j,k,l)
                        pre_vol(j,k,l)=post_vol(j,k,l)+vol_flux_y(j  ,k+1,l  )-vol_flux_y(j,k,l)

                        ENDDO
                ENDDO
            ENDDO
        ENDDO
    ENDDO
ENDDO
!$OMP END DO


 iter=(z_max-z_min+1+4)*(y_max-y_min+1+4)*(x_max-x_min+1+4)
 tmpFlop=iter*4
 tmpMem=iter*8
 


ELSEIF(sweep_number.EQ.2.AND..NOT.advect_x)THEN ! Z first
!$OMP DO collapse(3)
DO bl=z_min-2,nb_z,BLOCK_SIZE_z
    DO bk=y_min-2,nb_y,BLOCK_SIZE_y
      DO bj=x_min-2,nb_x,BLOCK_SIZE_x
        DO l=bl,min(nb_z,BLOCK_SIZE_z+bl-1)
            DO k=bk,min(nb_y,BLOCK_SIZE_y+bk-1)
                DO j=bj,min(nb_x,BLOCK_SIZE_x+bj-1)
                        post_vol(j,k,l)=volume(j,k,l) +vol_flux_x(j+1,k  ,l  )-vol_flux_x(j,k,l)
                        pre_vol(j,k,l)=post_vol(j,k,l)+vol_flux_y(j  ,k+1,l  )-vol_flux_y(j,k,l)

                        ENDDO
                ENDDO
            ENDDO
        ENDDO
    ENDDO
ENDDO
!$OMP END DO


  iter=(z_max-z_min+1+4)*(y_max-y_min+1+4)*(x_max-x_min+1+4)
  tmpFlop=iter*4
  tmpMem=iter*8

  
  ELSEIF(sweep_number.EQ.3.AND.direction.EQ.1)THEN ! z first
!$OMP DO collapse(3)
DO bl=z_min-2,nb_z,BLOCK_SIZE_z
    DO bk=y_min-2,nb_y,BLOCK_SIZE_y
      DO bj=x_min-2,nb_x,BLOCK_SIZE_x
        DO l=bl,min(nb_z,BLOCK_SIZE_z+bl-1)
            DO k=bk,min(nb_y,BLOCK_SIZE_y+bk-1)
                DO j=bj,min(nb_x,BLOCK_SIZE_x+bj-1)
                  post_vol(j,k,l)=volume(j,k,l)
                  pre_vol(j,k,l)=post_vol(j,k,l)+vol_flux_x(j+1,k  ,l  )-vol_flux_x(j,k,l)
                ENDDO
            ENDDO
        ENDDO
     ENDDO
   ENDDO
ENDDO
!$OMP END DO


  iter=(z_max-z_min+1+4)*(y_max-y_min+1+4)*(x_max-x_min+1+4)
  tmpFlop=iter*2
  tmpMem=iter*6


  ELSEIF(sweep_number.EQ.3.AND.direction.EQ.3)THEN ! x first
!$OMP DO collapse(3)
DO bl=z_min-2,nb_z,BLOCK_SIZE_z
    DO bk=y_min-2,nb_y,BLOCK_SIZE_y
      DO bj=x_min-2,nb_x,BLOCK_SIZE_x
        DO l=bl,min(nb_z,BLOCK_SIZE_z+bl-1)
            DO k=bk,min(nb_y,BLOCK_SIZE_y+bk-1)
                DO j=bj,min(nb_x,BLOCK_SIZE_x+bj-1)
                    post_vol(j,k,l)=volume(j,k,l)
                    pre_vol(j,k,l)=post_vol(j,k,l)+vol_flux_z(j  ,k  ,l+1)-vol_flux_z(j,k,l)
                ENDDO
                ENDDO
            ENDDO
        ENDDO
    ENDDO
ENDDO

!$OMP END DO



  iter=(z_max-z_min+1+4)*(y_max-y_min+1+4)*(x_max-x_min+1+4)
  tmpFlop=iter*2
  tmpMem=iter*6


  ENDIF
  

  flop=flop+tmpFlop
  mem=mem+tmpMem 
  tmpFlop=0  
  tmpMem=0
  


!----------------------------------- PRE and POST Volume calculation complete ----------------------------------

nb_x=x_max+2
nb_y=y_max+1
nb_z=z_max+1

  IF(direction.EQ.1)THEN
    IF(which_vel.EQ.1) THEN
!$OMP DO collapse(3)
    DO bl=z_min,nb_z,BLOCK_SIZE_z
        DO bk=y_min,nb_y,BLOCK_SIZE_y
          DO bj=x_min-2,nb_x,BLOCK_SIZE_x
            DO l=bl,min(nb_z,BLOCK_SIZE_z+bl-1)
                DO k=bk,min(nb_y,BLOCK_SIZE_y+bk-1)
                    DO j=bj,min(nb_x,BLOCK_SIZE_x+bj-1)
                        ! Find staggered mesh mass fluxes, nodal masses and volumes.
                            node_flux(j,k,l)=0.125_8*(mass_flux_x(j  ,k-1,l  )+mass_flux_x(j  ,k,l  )  &
                                     +mass_flux_x(j+1,k-1,l  )+mass_flux_x(j+1,k,l  )  &
                                     +mass_flux_x(j  ,k-1,l-1)+mass_flux_x(j  ,k,l-1)  &
                                     +mass_flux_x(j+1,k-1,l-1)+mass_flux_x(j+1,k,l-1))


                            !print *,"node_flux(,",j,",",k,",",l,"):",node_flux(j,k,l)
                        ENDDO
                   ENDDO
                ENDDO
            ENDDO
        ENDDO
      ENDDO
!$OMP END DO


  iter=(z_max-z_min+2)*(y_max-y_min+2)*(x_max-x_min+1+4)
  tmpFlop=iter*8
  tmpMem=iter*9
   flop=flop+tmpFlop !sc 27664
   mem=mem+tmpMem
   tmpFlop=0 
   tmpMem=0


!And do I need to calc the node mass for all 3 directions, or just once?

nb_x=x_max+2
nb_y=y_max+1
nb_z=z_max+1


!$OMP DO collapse(3)
DO bl=z_min,nb_z,BLOCK_SIZE_z
    DO bk=y_min,nb_y,BLOCK_SIZE_y
      DO bj=x_min,nb_x,BLOCK_SIZE_x
        DO l=bl,MIN(nb_z,BLOCK_SIZE_z+bl-1)
             DO k=bk,MIN(BLOCK_SIZE_y+bk-1,nb_y)
                DO j=bj,min(nb_x,BLOCK_SIZE_x+bj-1)
                ! Staggered cell mass post advection
                node_mass_post(j,k,l)=0.125_8*(density1(j  ,k-1,l  )*post_vol(j  ,k-1,l  )                   &
                                          +density1(j  ,k  ,l  )*post_vol(j  ,k  ,l  )                   &
                                          +density1(j-1,k-1,l  )*post_vol(j-1,k-1,l  )                   &
                                          +density1(j-1,k  ,l  )*post_vol(j-1,k  ,l  )                   &
                                          +density1(j  ,k-1,l-1)*post_vol(j  ,k-1,l-1)                   &
                                          +density1(j  ,k  ,l-1)*post_vol(j  ,k  ,l-1)                   &
                                          +density1(j-1,k-1,l-1)*post_vol(j-1,k-1,l-1)                   &
                                          +density1(j-1,k  ,l-1)*post_vol(j-1,k  ,l-1))


                !print *,"node_mass_post(,",j,",",k,",",l,"):",node_mass_post(j,k,l)
                ENDDO
                ENDDO
           ENDDO
       ENDDO
   ENDDO
ENDDO
!$OMP END DO



      iter=(z_max-z_min+2)*(y_max-y_min+2)*(x_max-x_min+1+3)
      tmpFlop=iter*16
      tmpMem=iter*17
      flop=flop+tmpFlop !sc 48464
      mem=mem+tmpMem
      tmpFlop=0
      tmpMem=0


nb_x=x_max+2
nb_y=y_max+1
nb_z=z_max+1


!$OMP DO collapse(3)
DO bl=z_min,nb_z,BLOCK_SIZE_z
    DO bk=y_min,nb_y,BLOCK_SIZE_y
        DO bj=x_min-1,nb_x,BLOCK_SIZE_x
            DO l=bl,min(nb_z,BLOCK_SIZE_z+bl-1)
                DO k=bk,min(nb_y,BLOCK_SIZE_y+bk-1)
                    DO j=bj,min(nb_x,BLOCK_SIZE_x+bj-1)
                    ! Staggered cell mass pre advection
                       node_mass_pre(j,k,l)=node_mass_post(j,k,l)-node_flux(j-1,k,l)+node_flux(j,k,l)

                       !print *,"node_mass_pre(,",j,",",k,",",l,"):",node_mass_pre(j,k,l)
                ENDDO
                ENDDO
            ENDDO
        ENDDO
    ENDDO
ENDDO
!$OMP END DO



     iter=(z_max-z_min+2)*(y_max-y_min+2)*(x_max-x_min+4)
     tmpFlop=iter*2
     tmpMem=iter*4
     flop=flop+tmpFlop !sc 51064
     mem=mem+tmpMem
     tmpFlop=0
     tmpMem=0



    ENDIF

!$OMP DO PRIVATE(upwind,downwind,donor,dif,sigma,width,limiter,vdiffuw,vdiffdw,auw,adw,wind) REDUCTION (+:limiter_count)
    DO l=z_min,z_max+1
      DO k=y_min,y_max+1
        DO j=x_min-1,x_max+1
          IF(node_flux(j,k,l).LT.0.0)THEN
            upwind=j+2
            donor=j+1
            downwind=j
            dif=donor
          ELSE
            upwind=j-1
            donor=j
            downwind=j+1
            dif=upwind
          ENDIF
          sigma=ABS(node_flux(j,k,l))/(node_mass_pre(donor,k,l))
          width=celldx(j)

          !----------- velocity difference, upwind and downwind -------------------------------------------------------------------

          vdiffuw=vel1(donor,k,l)-vel1(upwind,k,l)
          vdiffdw=vel1(downwind,k,l)-vel1(donor,k,l)

          !------------ start limiter calculation --------------------------------------------------------------------------------
          limiter=0.0
          IF(vdiffuw*vdiffdw.GT.0.0)THEN
            ! Sid - absolute value of the velocity difference, upwind and downwind
            auw=ABS(vdiffuw)
            adw=ABS(vdiffdw)

            ! Sid - set wind = 1.0
            wind=1.0_8
            IF(vdiffdw.LE.0.0) wind=-1.0_8

            !Sid - Wind is now set            
            !sigma is the flux/pre mass, tmp1= 2/sigma
            !Sid - absolute value of downwind divided by width
            !Sid - Sigma + 1
            !Sid - absolute value of upwind/celldx           
            limiter=wind*MIN(width*((2.0_8-sigma)*adw/width+(1.0_8+sigma)*auw/celldx(dif))/6.0_8,auw,adw)
            limiter_count=limiter_count+1
          ENDIF
          advec_vel_s=vel1(donor,k,l)+(1.0-sigma)*limiter
          mom_flux(j,k,l)=advec_vel_s*node_flux(j,k,l)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO


   iter=(z_max-z_min+2)*(y_max-y_min+2)*(x_max-x_min+3)
   tmpFlop=tmpFlop+(iter*10)+(limiter_count*19)
   tmpMem=tmpMem+(iter*24)+(limiter_count*10)
   flop=flop+tmpFlop
   mem=mem+tmpMem
   tmpFlop=0
   tmpMem=0


!------------------------------------ Limiter calculation complete ------------------------------------------------------
nb_x=x_max+1
nb_y=y_max+1
nb_z=z_max+1

!$OMP DO collapse(3)
DO bl=z_min,nb_z,BLOCK_SIZE_z
    DO bk=y_min,nb_y,BLOCK_SIZE_y
       DO bj=x_min,nb_x,BLOCK_SIZE_x
            DO l=bl,min(nb_z,BLOCK_SIZE_z+bl-1)
                DO k=bk,min(nb_y,BLOCK_SIZE_y+bk-1)
                   DO j=bj,min(nb_x,BLOCK_SIZE_x+bj-1)
                            vel1 (j,k,l)=(vel1 (j,k,l)*node_mass_pre(j,k,l)+mom_flux(j-1,k,l)-mom_flux(j,k,l))/node_mass_post(j,k,l)

                            !print *,"vel1(",j,",",k,",",l,"):",vel1(j,k,l)
                        ENDDO
                ENDDO
            ENDDO
        ENDDO
    ENDDO
ENDDO
!$OMP END DO



  iter=(z_max-z_min+2)*(y_max-y_min+2)*(x_max-x_min+2)
  tmpFlop=iter*7
  tmpMem=iter*6
  flop=flop+tmpFlop
  mem=mem+tmpMem
  tmpFlop=0
  tmpMem=0


!------------------------------ Sid - Velocity calculation complete ------------------------------------------------------
!------------------------------ repeat the same logic below for direction 2 and direction 3 ------------------------------

  ELSEIF(direction.EQ.2)THEN
    IF(which_vel.EQ.1)THEN


nb_x=x_max+1
nb_y=y_max+2
nb_z=z_max+1

!$OMP DO collapse(3)
DO bl=z_min,nb_z,BLOCK_SIZE_z
    DO bk=y_min-2,nb_y,BLOCK_SIZE_y

      DO bj=x_min,nb_x,BLOCK_SIZE_x
        DO l=bl,min(nb_z,BLOCK_SIZE_z+bl-1)
            DO k=bk,min(nb_y,BLOCK_SIZE_y+bk-1)
                DO j=bj,min(nb_x,BLOCK_SIZE_x+bj-1)
               ! Find staggered mesh mass fluxes and nodal masses and volumes.
                node_flux(j,k,l)=0.125_8*(mass_flux_y(j-1,k  ,l  )+mass_flux_y(j  ,k  ,l  ) &
                                     +mass_flux_y(j-1,k+1,l  )+mass_flux_y(j  ,k+1,l  ) &
                                     +mass_flux_y(j-1,k  ,l-1)+mass_flux_y(j  ,k  ,l-1) &
                                     +mass_flux_y(j-1,k+1,l-1)+mass_flux_y(j  ,k+1,l-1))

                ENDDO

                ENDDO
            ENDDO
        ENDDO
    ENDDO
ENDDO
!$OMP END DO


      iter=(z_max-z_min+2)*(y_max-y_min+5)*(x_max-x_min+2)
      tmpFlop=iter*8
      tmpMem=iter*9
      mem=mem+tmpMem
      flop=flop+tmpFlop
      tmpFlop=0
      tmpMem=0

nb_x=x_max+1
nb_y=y_max+2
nb_z=z_max+1


!$OMP DO collapse(3)
DO bl=z_min,nb_z,BLOCK_SIZE_z
    DO bk=y_min-1,nb_y,BLOCK_SIZE_y
      DO bj=x_min,nb_x,BLOCK_SIZE_x
        DO l=bl,min(nb_z,BLOCK_SIZE_z+bl-1)
            DO k=bk,min(nb_y,BLOCK_SIZE_y+bk-1)
                DO j=bj,min(nb_x,BLOCK_SIZE_x+bj-1)
                       node_mass_post(j,k,l)=0.125_8*(density1(j  ,k-1,l  )*post_vol(j  ,k-1,l  )                     &
                                          +density1(j  ,k  ,l  )*post_vol(j  ,k  ,l  )                     &
                                          +density1(j-1,k-1,l  )*post_vol(j-1,k-1,l  )                     &
                                          +density1(j-1,k  ,l  )*post_vol(j-1,k  ,l  )                     &
                                          +density1(j  ,k-1,l-1)*post_vol(j  ,k-1,l-1)                     &
                                          +density1(j  ,k  ,l-1)*post_vol(j  ,k  ,l-1)                     &
                                          +density1(j-1,k-1,l-1)*post_vol(j-1,k-1,l-1)                     &
                                          +density1(j-1,k  ,l-1)*post_vol(j-1,k  ,l-1))



                     ! print *,"node_mass_post(,",j,",",k,",",l,"):",node_mass_post(j,k,l)
                ENDDO
                ENDDO
            ENDDO
        ENDDO
    ENDDO
ENDDO
!$OMP END DO


      iter=(z_max-z_min+2)*(y_max-y_min+4)*(x_max-x_min+2)
      tmpFlop=iter*16
      tmpMem=iter*17
      flop=flop+tmpFlop
      mem=mem+tmpMem
      tmpMem=0
      tmpFlop=0


nb_x=x_max+1
nb_y=y_max+2
nb_z=z_max+1


!$OMP DO collapse(3)
DO bl=z_min,nb_z,BLOCK_SIZE_z
    DO bk=y_min-1,nb_y,BLOCK_SIZE_y
      DO bj=x_min,nb_x,BLOCK_SIZE_x
        DO l=bl,min(nb_z,BLOCK_SIZE_z+bl-1)
            DO k=bk,min(nb_y,BLOCK_SIZE_y+bk-1)
                DO j=bj,min(nb_x,BLOCK_SIZE_x+bj-1)
                    ! Staggered cell mass pre advection
                            node_mass_pre(j,k,l)=node_mass_post(j,k,l)-node_flux(j,k-1,l)+node_flux(j,k,l)
                ENDDO
                ENDDO
            ENDDO
        ENDDO
    ENDDO
ENDDO
!$OMP END DO


    iter=(z_max-z_min+2)*(y_max-y_min+4)*(x_max-x_min+2)
    tmpFlop=iter*2
    tmpMem=iter*4
    flop=flop+tmpFlop
    mem=mem+tmpMem
    tmpFlop=0
    tmpMem=0


ENDIF

limiter_count=0
!$OMP DO PRIVATE(upwind,donor,downwind,dif,sigma,width,limiter,vdiffuw,vdiffdw,auw,adw,wind) REDUCTION (+:limiter_count)
    DO l=z_min,z_max+1
      DO k=y_min-1,y_max+1
        DO j=x_min,x_max+1
          IF(node_flux(j,k,l).LT.0.0)THEN
            upwind=k+2
            donor=k+1
            downwind=k
            dif=donor
          ELSE
            upwind=k-1
            donor=k
            downwind=k+1
            dif=upwind
          ENDIF

          sigma=ABS(node_flux(j,k,l))/(node_mass_pre(j,donor,l))
          width=celldy(k)
          vdiffuw=vel1(j,donor,l)-vel1(j,upwind,l)
          vdiffdw=vel1(j,downwind,l)-vel1(j,donor,l)
          limiter=0.0

          IF(vdiffuw*vdiffdw.GT.0.0)THEN
            auw=ABS(vdiffuw)
            adw=ABS(vdiffdw)
            wind=1.0_8
            IF(vdiffdw.LE.0.0) wind=-1.0_8
            limiter=wind*MIN(width*((2.0_8-sigma)*adw/width+(1.0_8+sigma)*auw/celldy(dif))/6.0_8,auw,adw)
            limiter_count=limiter_count+1
          ENDIF
          advec_vel_s=vel1(j,donor,l)+(1.0_8-sigma)*limiter
          mom_flux(j,k,l)=advec_vel_s*node_flux(j,k,l)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO


  iter=(z_max-z_min+2)*(y_max-y_min+3)*(x_max-x_min+2)
  tmpFlop=tmpFlop+(iter*9)+(limiter_count*19)
  tmpMem=tmpMem+(iter*24)+(limiter_count*20)
  flop=flop+tmpFlop
  mem=mem+tmpMem
  tmpFlop=0
  tmpMem=0

 

!###################################

nb_x=x_max+1
nb_y=y_max+1
nb_z=z_max+1

!$OMP DO collapse(3)

DO bl=z_min,nb_z,BLOCK_SIZE_z
   DO bk=y_min,nb_y,BLOCK_SIZE_y
      DO bj=x_min,nb_x,BLOCK_SIZE_x
          DO l=bl,min(nb_z,BLOCK_SIZE_z+bl-1)
              DO k=bk,min(nb_y,BLOCK_SIZE_y+bk-1)
                  DO j=bj,min(nb_x,BLOCK_SIZE_x+bj-1)
                        vel1 (j,k,l)=(vel1(j,k,l)*node_mass_pre(j,k,l)+mom_flux(j,k-1,l)-mom_flux(j,k,l))/node_mass_post(j,k,l)
                  ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
  ENDDO
!$OMP END DO


  iter=(z_max-z_min+2)*(y_max-y_min+2)*(x_max-x_min+2)
  tmpFlop=iter*7
  tmpMem=iter*6
  flop=flop+tmpFlop
  mem=mem+tmpMem
  tmpFlop=0
  tmpMem=0


  ELSEIF(direction.EQ.3)THEN
    IF(which_vel.EQ.1) THEN




nb_x=x_max+1
nb_y=y_max+1
nb_z=z_max+2

!$OMP DO collapse(3)
DO bl=z_min-2,nb_z,BLOCK_SIZE_z
    DO bk=y_min,nb_y,BLOCK_SIZE_y
      DO bj=x_min,nb_x,BLOCK_SIZE_x
        DO l=bl,min(nb_z,BLOCK_SIZE_z+bl-1)
            DO k=bk,min(nb_y,BLOCK_SIZE_y+bk-1)
                DO j=bj,min(nb_x,BLOCK_SIZE_x+bj-1)
                    ! Find staggered mesh mass fluxes and nodal masses and volumes.
                    node_flux(j,k,l)=0.125_8*(mass_flux_z(j-1,k  ,l  )+mass_flux_z(j  ,k  ,l  ) &
                                     +mass_flux_z(j-1,k  ,l+1)+mass_flux_z(j  ,k  ,l+1) &
                                     +mass_flux_z(j-1,k-1,l  )+mass_flux_z(j  ,k-1,l  ) &
                                     +mass_flux_z(j-1,k-1,l+1)+mass_flux_z(j  ,k-1,l+1))
                ENDDO
                ENDDO
            ENDDO
        ENDDO
    ENDDO
ENDDO
!$OMP END DO


      iter=(z_max-z_min+5)*(y_max-y_min+2)*(x_max-x_min+2)
      tmpFlop=iter*8
      tmpMem=iter*9
      flop=flop+tmpFlop
      mem=mem+tmpMem
      tmpFlop=0
      tmpMem=0


nb_x=x_max+1
nb_y=y_max+1
nb_z=z_max+2


!$OMP DO collapse(3)
DO bl=z_min-1,nb_z,BLOCK_SIZE_z
    DO bk=y_min,nb_y,BLOCK_SIZE_y
      DO bj=x_min,nb_x,BLOCK_SIZE_x
        DO l=bl,min(nb_z,BLOCK_SIZE_z+bl-1)
            DO k=bk,min(nb_y,BLOCK_SIZE_y+bk-1)
                DO j=bj,min(nb_x,BLOCK_SIZE_x+bj-1)
                    node_mass_post(j,k,l)=0.125_8*(density1(j  ,k-1,l  )*post_vol(j  ,k-1,l  )                     &
                                          +density1(j  ,k  ,l  )*post_vol(j  ,k  ,l  )                     &
                                          +density1(j-1,k-1,l  )*post_vol(j-1,k-1,l  )                     &
                                          +density1(j-1,k  ,l  )*post_vol(j-1,k  ,l  )                     &
                                          +density1(j  ,k-1,l-1)*post_vol(j  ,k-1,l-1)                     &
                                          +density1(j  ,k  ,l-1)*post_vol(j  ,k  ,l-1)                     &
                                          +density1(j-1,k-1,l-1)*post_vol(j-1,k-1,l-1)                     &
                                          +density1(j-1,k  ,l-1)*post_vol(j-1,k  ,l-1))
                ENDDO
                ENDDO
            ENDDO
        ENDDO
    ENDDO
ENDDO
!$OMP END DO


      iter=(z_max-z_min+4)*(y_max-y_min+2)*(x_max-x_min+2)
      tmpFlop=iter*16
      tmpMem=iter*17
      mem=mem+tmpMem
      flop=flop+tmpFlop
      tmpFlop=0
      tmpMem=0

nb_x=x_max+1
nb_y=y_max+1
nb_z=z_max+2


!$OMP DO collapse(3)
DO bl=z_min-1,nb_z,BLOCK_SIZE_z
    DO bk=y_min,nb_y,BLOCK_SIZE_y
      DO bj=x_min,nb_x,BLOCK_SIZE_x
        DO l=bl,min(nb_z,BLOCK_SIZE_z+bl-1)
            DO k=bk,min(nb_y,BLOCK_SIZE_y+bk-1)
                DO j=bj,min(nb_x,BLOCK_SIZE_x+bj-1)
                        ! Staggered cell mass pre advection
                           node_mass_pre(j,k,l)=node_mass_post(j,k,l)-node_flux(j,k,l-1)+node_flux(j,k,l)
                ENDDO
                ENDDO
            ENDDO
        ENDDO
    ENDDO
ENDDO
!$OMP END DO


    iter=(z_max-z_min+4)*(y_max-y_min+2)*(x_max-x_min+2)
    tmpFlop=iter*2
    tmpMem=iter*4
    flop=flop+tmpFlop
    mem=mem+tmpMem
    tmpFlop=0
    tmpMem=0


ENDIF

limiter_count=0
!$OMP DO PRIVATE(upwind,donor,downwind,dif,sigma,width,limiter,vdiffuw,vdiffdw,auw,adw,wind) REDUCTION (+:limiter_count)
    DO l=z_min-1,z_max+1
      DO k=y_min,y_max+1
        DO j=x_min,x_max+1
          IF(node_flux(j,k,l).LT.0.0)THEN
            upwind=l+2
            donor=l+1
            downwind=l
            dif=donor
          ELSE
            upwind=l-1
            donor=l
            downwind=l+1
            dif=upwind
          ENDIF

          sigma=ABS(node_flux(j,k,l))/(node_mass_pre(j,k,donor))
          width=celldz(l)
          vdiffuw=vel1(j,k,donor)-vel1(j,k,upwind)
          vdiffdw=vel1(j,k,downwind)-vel1(j,k,donor)
          limiter=0.0
          IF(vdiffuw*vdiffdw.GT.0.0)THEN
            auw=ABS(vdiffuw)
            adw=ABS(vdiffdw)
            wind=1.0_8
            IF(vdiffdw.LE.0.0) wind=-1.0_8
            limiter=wind*MIN(width*((2.0_8-sigma)*adw/width+(1.0_8+sigma)*auw/celldz(dif))/6.0_8,auw,adw)
            limiter_count=limiter_count+1
          ENDIF
          advec_vel_s=vel1(j,k,donor)+(1.0_8-sigma)*limiter
          mom_flux(j,k,l)=advec_vel_s*node_flux(j,k,l)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO



  iter=(z_max-z_min+2)*(y_max-y_min+2)*(x_max-x_min+3)
  tmpFlop=tmpFlop+(iter*10)+(limiter_count*18)
  tmpMem=tmpMem+(iter*24)+(limiter_count*20)
  flop=flop+tmpFlop
  mem=mem+tmpMem
  tmpFlop=0
  tmpMem=0


!###############################################

nb_x=x_max+1
nb_y=y_max+1
nb_z=z_max+1

!$OMP DO collapse(3)

DO bl=z_min,nb_z,BLOCK_SIZE_z
    DO bk=y_min,nb_y,BLOCK_SIZE_y
      DO bj=x_min,nb_x,BLOCK_SIZE_x
        DO l=bl,min(nb_z,BLOCK_SIZE_z+bl-1)
            DO k=bk,min(nb_y,BLOCK_SIZE_y+bk-1)
                DO j=bj,min(nb_x,BLOCK_SIZE_x+bj-1)
                        vel1 (j,k,l)=(vel1(j,k,l)*node_mass_pre(j,k,l)+mom_flux(j,k,l-1)-mom_flux(j,k,l))/node_mass_post(j,k,l)
                ENDDO
                ENDDO
            ENDDO
        ENDDO
    ENDDO
ENDDO
!$OMP END DO
!###############################################


  iter=(z_max-z_min+2)*(y_max-y_min+2)*(x_max-x_min+2)
  tmpFlop=iter*7
  tmpMem=iter*6


  ENDIF


  flop=flop+tmpFlop
  mem=mem+tmpMem
  
!$OMP END PARALLEL
END SUBROUTINE advec_mom_kernel

END MODULE advec_mom_kernel_mod
