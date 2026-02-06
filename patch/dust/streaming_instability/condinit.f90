!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use amr_commons
  use hydro_parameters
  use hydro_commons
  use cooling_module      , only : kb,mh
  use radiation_parameters,only:mu_gas
  use cloud_module
  use units_commons
  
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar+3)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:4): d.u,d.v,d.w, U(i,5): E, U(i,6:8): Bleft, 
  ! U(i,nvar+1:nvar+3): Bright
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:4):u,v,w, Q(i,5): P, Q(i,6:8): Bleft, 
  ! Q(i,nvar+1:nvar+3): Bright
  ! If nvar > 8, remaining variables (9:nvar) are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer::ivar,i,j,k,idust,igroup
  real(dp),dimension(1:nvector,1:nvar+3),save::q   ! Primitive variables
  real(dp)::pi,xx,yy,xn,yn,zn
  real(dp)::mag_norm,scale_dcol

  real(dp),save:: first
  real(dp),dimension(1:3,1:100,1:100,1:100),save::q_idl
  real(dp),save::vx_tot,vy_tot,vz_tot,vx2_tot,vy2_tot,vz2_tot,vx,vy,vz,v_rms
  integer,save:: n_size
  integer:: ind_i, ind_j, ind_k
  real(dp),save:: ind,seed1,seed2,seed3,xi,yi,zi
  real(dp):: n_total
  real(dp):: u_rf,tempe,cs,zcoro,k_corona,coro_fact,cs_k,dvx,drho
  real(dp):: Height0,height0dust,rhocen,Temp0,rhosurf,sigma_dust,radiation_source
  
  real(dp),dimension(1:ndim) :: vtur
#if NDUST>0
  real(dp),dimension(1:ndust):: dustMRN
#endif
  real(dp):: epsilon_0,sum_dust
#if NDUST>0
  epsilon_0 = dust_ratio(1)
#endif

  pi=ACOS(-1.0d0)


  scale_dcol= scale_m/scale_l**2.0


  cs = hover_r*omega_shear*radii0


  do i=1,nn

     xn = x(i,1) !- 0.5*boxlen
     yn = x(i,2) !- 0.5*boxlen
     zn = x(i,3) !- 0.5*boxlen
      
#if NDUST>0
     if(mrn) call init_dust_ratio(epsilon_0, dustMRN)
#endif     
     !Bx component 
     q(i,6     ) = 0.0d0
     q(i,nvar+1) =0.0d0

     !By component
     q(i,7     ) = 0.0d0
     q(i,nvar+2) = 0.0d0

     !Bz component
     q(i,8     ) =  0.0d0
     q(i,nvar+3) =  0.0d0

     !densite
     ky_0=ky_0!/((hover_r*radii0)**2/radii0)
     kz_0=kz_0!/((hover_r*radii0)**2/radii0)
     call get_vxturb(turb_perc,rho0,drho)

     q(i,1) = rho0*(1.0d0+amplitude_pert*(cos(2.0*pi*(ky_0*abs(yn)/boxlen+kz_0*abs(zn)/boxlen)))) +drho

#if NDUST>0
     sum_dust=0.0d0
     do idust =1,ndust

        q(i, firstindex_ndust+idust)=rhod0
        !sum_dust = sum_dust +drho
        sum_dust = sum_dust + q(i, firstindex_ndust+idust)
                
        end do   
#endif
        q(i,5) = q(i,1)*cs**2.0
        q(i,1)=q(i,1)+ sum_dust
    
     call get_vxturb(0.0d0,cs,dvx)
     q(i,2) = q_shear*Omega_shear*boxlen*((x(i,2) - 0.5*boxlen)/boxlen) + dvx
     call get_vxturb(0.0d0,cs,dvx)
     q(i,3) = dvx
     call get_vxturb(0.0d0,cs,dvx)
     q(i,4) =dvx
     

  end do
  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
  ! kinetic energy
  u(1:nn,5)=0.0d0
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,2)**2
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,3)**2
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,4)**2
  ! pressure -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+q(1:nn,5)/(gamma-1.0d0)
  ! magnetic energy -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,6)+q(1:nn,nvar+1))**2
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,7)+q(1:nn,nvar+2))**2
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,8)+q(1:nn,nvar+3))**2
  u(1:nn,6:8)=q(1:nn,6:8)
  u(1:nn,nvar+1:nvar+3)=q(1:nn,nvar+1:nvar+3)
  u(1:nn,nvar)=q(1:nn,5)/(gamma-1.0d0)

  ! passive scalars
  do ivar=9+nener,nvar
     u(1:nn,ivar)= q(1:nn,ivar)
  end do

end subroutine condinit
!================================================================
!================================================================
!================================================================
!================================================================
subroutine velana(x,v,dx,t,ncell)
  use amr_parameters
  use hydro_parameters  
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp)::t                             ! Current time
  real(dp),dimension(1:nvector,1:3)::v    ! Velocity field
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the user defined velocity fields.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! v(i,1:3) is the imposed 3-velocity in user units.
  !================================================================
  integer::i
  real(dp)::xx,yy,zz,vx,vy,vz,rr,tt,omega,aa,twopi

  ! Add here, if you wish, some user-defined initial conditions
  aa=1.0
  twopi=2d0*ACOS(-1d0)
  do i=1,ncell

     xx=x(i,1)
#if NDIM > 1
     yy=x(i,2)
#endif
#if NDIM > 2
     zz=x(i,3)
#endif
     ! ABC
     vx=aa*(cos(twopi*yy)+sin(twopi*zz))
     vy=aa*(sin(twopi*xx)+cos(twopi*zz))
     vz=aa*(cos(twopi*xx)+sin(twopi*yy))

     v(i,1)=vx
#if NDIM > 1
     v(i,2)=vy
#endif
#if NDIM > 2
     v(i,3)=vz
#endif
  end do


end subroutine velana
subroutine get_vturb(vrms,cs,v_turb)

  use amr_commons, only:myid,ncpu
  use pm_commons, only:localseed,IRandNumSize,iseed
  use random

  implicit none

  integer :: i
  integer ,dimension(1:ncpu,1:IRandNumSize)    :: allseed
  double precision ::  vrms, cs, vel
  double precision, dimension(3) :: v_turb

  double precision :: u1, v1, u2, v2, theta, phi, x, y, z

#ifdef DEBUGRANDOM
  logical, save :: first=.true.
#endif

  if (localseed(1)==-1) then
     call rans(ncpu,iseed,allseed)
     localseed = allseed(myid,1:IRandNumSize)
  end if

  ! magnitude --> Gressel is v_rms = 0.01 * cs
  if (vrms .lt. 0.d0) then
    call gaussdev(localseed,vel)
    vel = vel * abs(vrms)*cs
  else
    call gaussdev(localseed,vel)
    vel = vel * vrms
  endif

  call ranf(localseed,v_turb(1))
  call ranf(localseed,v_turb(2))
  call ranf(localseed,v_turb(3))

  v_turb = (v_turb - 0.5d0)*2.d0
  v_turb = v_turb/sqrt(sum((v_turb(1:3))**2)) * vel

  ! NEED TO HAVE Cs OR Vrms in code units.

#ifdef DEBUGRANDOM
  if (myid == 1 .and. first) then
    open(42,file="test_gauss.dat",status='unknown')
    do i=1,10000
      call gaussdev(localseed,vel)
      write(42,*) vel
    enddo
    close(42)
    first = .false.
  endif
#endif

end subroutine get_vturb


subroutine get_vxturb(pert,vx,delvx)

  use amr_commons, only:myid,ncpu
  use pm_commons, only:localseed,IRandNumSize,iseed
  use random

  implicit none

  integer :: i
  integer ,dimension(1:ncpu,1:IRandNumSize)    :: allseed
  double precision ::  pert,vx,delvx,randno


#ifdef DEBUGRANDOM
  logical, save :: first=.true.
#endif

  if (localseed(1)==-1) then
     call rans(ncpu,iseed,allseed)
     localseed = allseed(myid,1:IRandNumSize)
  end if

  call ranf(localseed,randno)

  delvx = randno*pert* vx


end subroutine get_vxturb

