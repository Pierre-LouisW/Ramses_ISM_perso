!############################################################
!############################################################
!############################################################
!############################################################
subroutine boundana(x,u,vdu,dx,ibound,ncell,ilevel)
  use amr_parameters!, ONLY: dp,ndim,nvector
  use hydro_parameters!, ONLY: nvar,boundary_var
  use amr_commons
  use hydro_commons
  use cloud_module
  use cooling_module      , only : kb,mh
  use radiation_parameters,only:mu_gas
  use units_commons
  implicit none
  integer ::ibound                        ! Index of boundary region
  integer ::ncell                         ! Number of active cells
  real(dp)::dx,scale_dcol,temp0 ,enint,emag,ekin              ! Cell size
  real(dp),dimension(1:nvector,1:nvar+3)::u ! Conservative variables
    real(dp),dimension(1:nvector,1:ndust,1:ndim)::vdu

  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  integer,dimension(1:nvector)::cell_index,cell_index2
  integer,dimension(1:nvector)::cell_levl,cell_levl2
  real(dp),dimension(1:nvector,1:ndim)::y
  real(dp),dimension(1:nvector,1:ndim)::xtcell1,xtcell2,xtcell3
  !================================================================
  ! This routine generates boundary conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! U is in user units.
  ! ibound is the index of the boundary region defined in the namelist.
  !================================================================
  integer::ivar,i,idi,ilevel,idust

  real(dp)::dx_loc,diff,w1,w2,height0,cs,rhocen,zcoro,k_corona,pi
  dx_loc = 0.5d0**ilevel * boxlen
  dx = 0.5d0**ilevel
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! CC 01/12/2016
  ! y sont les coordonnees de la cellule du domaine correspondant a la cellule
  ! fantome, obtenues par periodicite et par decalage vers la gauche (boundary
  ! #2) ou la droite (boundary #1) du au shear
  ! Periodicite et shear sur la boundary 3 (en bas)
    pi=ACOS(-1.0d0)


  scale_dcol= scale_m/scale_l**2.0

  if(ibound.eq.1) then
     do i=1,ncell
        y(i,1)=(1.0D0/boxlen)*mod(x(i,1)+ q_shear*Omega_shear*boxlen*t,boxlen)
        if (y(i,1).lt.0) then
           y(i,1)=y(i,1)+1.0D0
        end if
        do idi=2,ndim
           y(i,idi)=(1.0D0/boxlen)*mod(x(i,idi),boxlen)
           if (y(i,idi).lt.0) then
              y(i,idi)=y(i,idi)+1.0D0
           end if
        end do
     end do

     call get_cell_index4(cell_index,cell_levl,y,xtcell1,ilevel,ncell)
    
     ! interpolation : no amr assumed here
     do i=1,ncell
        diff=y(i,1)-xtcell1(i,1)
        if(abs(diff)<dx/1000.) then
            diff=0d0
        else
            diff=diff/abs(diff)
        end if
        xtcell2(i,1)=mod(xtcell1(i,1)+dx*diff,1.d0)
        if(xtcell2(i,1)<0.) xtcell2(i,1)=xtcell2(i,1)+1.0D0
        xtcell2(i,2)=xtcell1(i,2)
        xtcell2(i,3)=xtcell1(i,3)
     end do

     call get_cell_index4(cell_index2,cell_levl2,xtcell2,xtcell3,ilevel,ncell)

     do i=1,ncell
        if(diff.eq.0d0) then
           w1 = 1.d0
        else if(abs(xtcell2(i,1)-y(i,1)) .le. dx) then
           w1 =  abs(xtcell2(i,1)-y(i,1))  / dx
        else
           w1 =  ( 1.-abs(xtcell2(i,1)-y(i,1))) / dx
        end if
     
        w2=1.-w1
        u(i,1) =  max(w1*uold(cell_index(i),1) + w2*uold(cell_index2(i),1) ,smallr)
        u(i,2) = w1*uold(cell_index(i),2) + w2*uold(cell_index2(i),2) - u(i,1)*q_shear*Omega_shear*boxlen
        vdu(i,:,1)=v_dust(cell_index(i),:,1)
        vdu(i,:,2)=v_dust(cell_index(i),:,2)
        vdu(i,:,3)=v_dust(cell_index(i),:,3)
        
        do ivar=3,nvar+3
           u(i,ivar)= w1*uold(cell_index(i),ivar) + w2*uold(cell_index2(i),ivar)
        end do
!!$        do ivar=firstindex_pscal,lastindex_pscal
!!$           u(i,ivar)=(w1*uold(cell_index(i),ivar)/uold(cell_index(i),1) + w2*uold(cell_index2(i),ivar)/uold(cell_index2(i),1))*u(i,1)
!!$        end do
!!$        do idust=firstindex_ndust+1,firstindex_ndust+ndust
!!$           u(i,idust)=u(i,idust)/uold(cell_index(i),1)*u(i,1)
!!$        end do
        ! Correct kinetic energy
        u(i,ndim+2) = u(i,ndim+2) - w1*0.5d0*( (uold(cell_index(i),2))**2 + (uold(cell_index(i),3))**2 + (uold(cell_index(i),4))**2 )/uold(cell_index(i),1)
        u(i,ndim+2) = u(i,ndim+2) - w2*0.5d0*( (uold(cell_index2(i),2))**2 + (uold(cell_index2(i),3))**2 + (uold(cell_index2(i),4))**2 )/uold(cell_index2(i),1)
        u(i,ndim+2) = u(i,ndim+2) + 0.5d0*(u(i,2)**2+u(i,3)**2+u(i,4)**2) / u(i,1)

        ! Correct magnetic energy
        u(i,ndim+2) = u(i,ndim+2) - w1*0.125D0*( (uold(cell_index(i),6)+uold(cell_index(i),nvar+1))**2 + (uold(cell_index(i),7)+uold(cell_index(i),nvar+2))**2 + (uold(cell_index(i),8)+uold(cell_index(i),nvar+3))**2 )
        u(i,ndim+2) = u(i,ndim+2) - w2*0.125D0*( (uold(cell_index2(i),6)+uold(cell_index2(i),nvar+1))**2 + (uold(cell_index2(i),7)+uold(cell_index2(i),nvar+2))**2 + (uold(cell_index2(i),8)+uold(cell_index2(i),nvar+3))**2 )
        u(i,ndim+2) = u(i,ndim+2) + 0.125D0*((u(i,6)+u(i,nvar+1))**2+(u(i,7)+u(i,nvar+2))**2+(u(i,8)+u(i,nvar+3))**2)

     end do 


     !now we must correct for div B. For that purpose we first grab the cells 
     !located in 0.5d0^levelmin and we propagate its By field
     do i=1,ncell
        y(i,2)=0.5d0*(0.5d0)**ilevel
        y(i,1)=(1.0D0/boxlen)*mod(x(i,1),boxlen)
        y(i,3)=(1.0D0/boxlen)*mod(x(i,3),boxlen)

        if (y(i,1).lt.0) then
           y(i,1)=y(i,1)+1.0D0
        end if
        if (y(i,3).lt.0) then
           y(i,3)=y(i,3)+1.0D0
        end if

     end do

     call get_cell_index4(cell_index,cell_levl,y,xtcell1,ilevel,ncell)

     do i=1,ncell
        !calculate the magnetic energy
        emag = 0.125D0*((u(i,6)+u(i,nvar+1))**2+(u(i,7)+u(i,nvar+2))**2+(u(i,8)+u(i,nvar+3))**2)
        !remove it from total energy
        u(i,ndim+2)=u(i,ndim+2)-emag

        !impose the continuity of Bperp
        u(i,nvar+2)=uold(cell_index(i),7)

        !u(i,6)=uold(cell_index(i),6)
        !u(i,8)=uold(cell_index(i),8)
        !u(i,nvar+1)=uold(cell_index(i),nvar+1)
        !u(i,nvar+3)=uold(cell_index(i),nvar+3)

        !now impose div B within the cell
        u(i,7)=u(i,nvar+2)-u(i,6)+u(i,nvar+1)-u(i,8)+u(i,nvar+3)
     end do 


     !now depending whether the ghost cell is located at -0.5d0^levelmin / 2 or at -0.5d0^levelmin / 2 * 3 
     !(in this latter case) we must grab the neigbour of the one located at -0.5d0^levelmin / 2 to get div B =0
     do i=1,ncell

        y(i,1)=(1.0D0/boxlen)*mod(x(i,1)+q_shear*Omega_shear*boxlen*t,boxlen)
        if (y(i,1).lt.0) then
           y(i,1)=y(i,1)+1.0D0
        end if

        !shift the cells from dx_loc
        y(i,2)=(1.0D0/boxlen)*mod(x(i,2)+dx_loc,boxlen)
        if (y(i,2).lt.0) then
           y(i,2)=y(i,2)+1.0D0
        end if

        y(i,3)=(1.0D0/boxlen)*mod(x(i,3),boxlen)
        if (y(i,3).lt.0) then
           y(i,3)=y(i,3)+1.0D0
        end if

     end do

     call get_cell_index4(cell_index,cell_levl,y,xtcell1,ilevel,ncell)

     do i=1,ncell
        !check whether the correction to the normal field  must be applied
        if(x(i,2) .lt. -dx_loc) then
           !impose the continuity of Bperp by adding the difference of the tangential field
!           u(i,nvar+2)=u(i,nvar+2)-uold(cell_index(i),6)+uold(cell_index(i),nvar+1)-uold(cell_index(i),8)+uold(cell_index(i),nvar+3)
           !apply the correction to the other component
!           u(i,7)=u(i,7)-uold(cell_index(i),6)+uold(cell_index(i),nvar+1)-uold(cell_index(i),8)+uold(cell_index(i),nvar+3)

           u(i,nvar+2)=u(i,nvar+2)-u(i,6)+u(i,nvar+1)-u(i,8)+u(i,nvar+3)
           u(i,7)=u(i,7)-u(i,6)+u(i,nvar+1)-u(i,8)+u(i,nvar+3)
        endif
        !calculate the magnetic energy
        emag = 0.125D0*((u(i,6)+u(i,nvar+1))**2+(u(i,7)+u(i,nvar+2))**2+(u(i,8)+u(i,nvar+3))**2)
        !add it back
        u(i,ndim+2)=u(i,ndim+2)+emag
     end do 


  end if




  ! Periodicite et shear sur la boundary 4 (en haut)
!  if(ibound.eq.4) then
  if(ibound.eq.2) then
     do i=1,ncell
        y(i,1)=(1.0D0/boxlen)*mod(x(i,1)-q_shear*Omega_shear*boxlen*t,boxlen)
        if (y(i,1).lt.0) then
           y(i,1)=y(i,1)+1.0D0
        end if
        do idi=2,ndim
           y(i,idi)=(1.0D0/boxlen)*mod(x(i,idi),boxlen)
           if (y(i,idi).lt.0) then
              y(i,idi)=y(i,idi)+1.0D0
           end if
        end do
     end do

     call get_cell_index4(cell_index,cell_levl,y,xtcell1,ilevel,ncell)

     ! interpolation : no amr assumed here
     do i=1,ncell
        diff=y(i,1)-xtcell1(i,1)
        if(abs(diff)<dx/1000.) then
            diff = 0d0
        else
            diff=diff/abs(diff) 
        end if
        xtcell2(i,1)=mod(xtcell1(i,1)+dx*diff,1.d0)
        if(xtcell2(i,1)<0.) xtcell2(i,1)=xtcell2(i,1)+1.0D0
        xtcell2(i,2)=xtcell1(i,2)
        xtcell2(i,3)=xtcell1(i,3)
     end do

     call get_cell_index4(cell_index2,cell_levl2,xtcell2,xtcell3,ilevel,ncell)

     do i=1,ncell
        if(diff.eq.0d0) then
           w1=1.d0
        else if(abs(xtcell2(i,1)-y(i,1)) .le. dx) then
           w1 =  abs(xtcell2(i,1)-y(i,1))  /  dx
        else
           w1 =  ( 1.-abs(xtcell2(i,1)-y(i,1)) ) /  dx
        end if
     
        w2=1.-w1
        
        u(i,1) = max(w1*uold(cell_index(i),1) + w2*uold(cell_index2(i),1) ,smallr)
        u(i,2) = w1*uold(cell_index(i),2) + w2*uold(cell_index2(i),2) + u(i,1)*q_shear*Omega_shear*boxlen
        vdu(i,:,1)=v_dust(cell_index(i),:,1)
        vdu(i,:,2)=v_dust(cell_index(i),:,2)
        vdu(i,:,3)=v_dust(cell_index(i),:,3)
        do ivar=3,nvar+3
           u(i,ivar)= w1*uold(cell_index(i),ivar) + w2*uold(cell_index2(i),ivar)
        end do
        !do ivar=firstindex_pscal,lastindex_pscal
        !   u(i,ivar)=(w1*uold(cell_index(i),ivar)/max(uold(cell_index(i),1),smallr) + w2*uold(cell_index2(i),ivar)/max(uold(cell_index2(i),1),smallr))*u(i,1)
        !end do
!!$        do ivar=firstindex_pscal,lastindex_pscal
!!$           u(i,ivar)=(w1*uold(cell_index(i),ivar)/uold(cell_index(i),1) + w2*uold(cell_index2(i),ivar)/uold(cell_index2(i),1))*u(i,1)
!!$        end do
!!$        do idust=firstindex_ndust+1,firstindex_ndust+ndust
!!$           u(i,idust)=u(i,idust)/uold(cell_index(i),1)*u(i,1)
!!$
!!$        end do
        ! Correct kinetic energy
        u(i,ndim+2) = u(i,ndim+2) - w1*0.5d0*( (uold(cell_index(i),2))**2 + (uold(cell_index(i),3))**2 + (uold(cell_index(i),4))**2 )/uold(cell_index(i),1)
        u(i,ndim+2) = u(i,ndim+2) - w2*0.5d0*( (uold(cell_index2(i),2))**2 + (uold(cell_index2(i),3))**2 + (uold(cell_index2(i),4))**2 )/uold(cell_index2(i),1)
        u(i,ndim+2) = u(i,ndim+2) + 0.5d0*(u(i,2)**2+u(i,3)**2+u(i,4)**2) / u(i,1)

        ! Correct magnetic energy
        u(i,ndim+2) = u(i,ndim+2) - w1*0.125D0*( (uold(cell_index(i),6)+uold(cell_index(i),nvar+1))**2 + (uold(cell_index(i),7)+uold(cell_index(i),nvar+2))**2 + (uold(cell_index(i),8)+uold(cell_index(i),nvar+3))**2 )
        u(i,ndim+2) = u(i,ndim+2) - w2*0.125D0*( (uold(cell_index2(i),6)+uold(cell_index2(i),nvar+1))**2 + (uold(cell_index2(i),7)+uold(cell_index2(i),nvar+2))**2 + (uold(cell_index2(i),8)+uold(cell_index2(i),nvar+3))**2 )
        u(i,ndim+2) = u(i,ndim+2) + 0.125D0*((u(i,6)+u(i,nvar+1))**2+(u(i,7)+u(i,nvar+2))**2+(u(i,8)+u(i,nvar+3))**2)

     end do 


     !now we must correct for div B. For that purpose we first grap the cells 
     !located in 0.5d0^levelmin and we propagate its By field
     do i=1,ncell
        y(i,2)=1.-0.5d0*(0.5d0)**ilevel
        y(i,1)=(1.0D0/boxlen)*mod(x(i,1),boxlen)
        y(i,3)=(1.0D0/boxlen)*mod(x(i,3),boxlen)

        if (y(i,1).lt.0) then
           y(i,1)=y(i,1)+1.0D0
        end if
        if (y(i,3).lt.0) then
           y(i,3)=y(i,3)+1.0D0
        end if

     end do

     call get_cell_index4(cell_index,cell_levl,y,xtcell1,ilevel,ncell)

     do i=1,ncell
        !calculate the magnetic energy
        emag = 0.125D0*((u(i,6)+u(i,nvar+1))**2+(u(i,7)+u(i,nvar+2))**2+(u(i,8)+u(i,nvar+3))**2)
        !remove it from total energy
        u(i,ndim+2)=u(i,ndim+2)-emag

        !impose the continuity of Bperp
        u(i,7)=uold(cell_index(i),nvar+2)
        !now impose div B within the cell
        u(i,nvar+2)=u(i,7)+u(i,6)-u(i,nvar+1)+u(i,8)-u(i,nvar+3)
     end do 


     !now depending whether the ghost cell is located at -0.5d0^levelmin / 2 or at -0.5d0^levelmin / 2 * 3 
     !(in this latter case) we must grasp the neigbour of the one located at -0.5d0^levelmin / 2 to get div B =0
     do i=1,ncell

        y(i,1)=(1.0D0/boxlen)*mod(x(i,1)+q_shear*Omega_shear*boxlen*t,boxlen)
        if (y(i,1).lt.0) then
           y(i,1)=y(i,1)+1.0D0
        end if

        !shift the cells from dx_loc
        y(i,2)=(1.0D0/boxlen)*mod(x(i,2)-dx_loc,boxlen)
        if (y(i,2).lt.0) then
           y(i,2)=y(i,2)+1.0D0
        end if

        y(i,3)=(1.0D0/boxlen)*mod(x(i,3),boxlen)
        if (y(i,3).lt.0) then
           y(i,3)=y(i,3)+1.0D0
        end if

     end do

     call get_cell_index4(cell_index,cell_levl,y,xtcell1,ilevel,ncell)

     do i=1,ncell
        !check whether the correction to the normal field  must be applied
        if(x(i,2) .gt. boxlen+dx_loc) then
           !impose the continuity of Bperp by adding the difference of the tangential field
           u(i,7)=u(i,7)+uold(cell_index(i),6)-uold(cell_index(i),nvar+1)+uold(cell_index(i),8)-uold(cell_index(i),nvar+3)
           !apply the correction to the other component
           u(i,nvar+2)=u(i,nvar+2)+uold(cell_index(i),6)-uold(cell_index(i),nvar+1)+uold(cell_index(i),8)-uold(cell_index(i),nvar+3)

        !u(i,7)=u(i,7)+u(i,6)-u(i,nvar+1)+u(i,8)-u(i,nvar+3)
        !u(i,nvar+2)=u(i,nvar+2)+u(i,6)-u(i,nvar+1)+u(i,8)-u(i,nvar+3)

        endif
        !calculate the magnetic energy
        emag = 0.125D0*((u(i,6)+u(i,nvar+1))**2+(u(i,7)+u(i,nvar+2))**2+(u(i,8)+u(i,nvar+3))**2)
        !add it back
        u(i,ndim+2)=u(i,ndim+2)+emag
     end do 

  end if



end subroutine boundana

subroutine get_cell_index4(cell_index,cell_levl,xpart,xtcell,ilevel,np)
  use amr_commons
  implicit none
  integer                                :: np,ilevel
  integer,dimension(1:nvector)           :: cell_index,cell_levl
  real(dp),dimension(1:nvector,1:3)      :: xpart
  ! This function returns the index of the cell, at maximum level
  ! ilevel, in which the input particle sits
  real(dp)                               :: xx,yy,zz
  integer                                :: i,j,ii,jj,kk,ind,iskip,igrid,ind_cell,igrid0,igrid_old

  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim)::xtcell

  integer::idim
  integer::ix,iy,iz
  real(dp)::dx
  real(dp),dimension(1:3)::skip_loc



  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)

  
  !if ((nx.eq.1).and.(ny.eq.1).and.(nz.eq.1)) then
  !else if ((nx.eq.3).and.(ny.eq.3).and.(nz.eq.3)) then
  !else
     !write(*,*)"nx=ny=nz != 1,3 is not supported."
     !call clean_stop
  !end if
  
  igrid0=son(1+icoarse_min+jcoarse_min*nx+kcoarse_min*nx*ny)
  do i=1,np
     xx = xpart(i,1)
     yy = xpart(i,2)
     zz = xpart(i,3)
     if( ((xx .le. 0) .or. (xx .ge. 1.)) .or. ((yy .le. 0) .or. (yy .ge. 1.)) .or. ((zz .le. 0) .or. (zz .ge. 1.)) ) then 
        cell_index(i)=-1.
     else 
        xx = xx + (nx-1.)/2.
        yy = yy + (ny-1.)/2.
        zz = zz + (nz-1.)/2.
        igrid=igrid0
        do j=1,ilevel
           ii=1; jj=1; kk=1
           if(xx<xg(igrid,1))ii=0
           if(yy<xg(igrid,2))jj=0
           if(zz<xg(igrid,3))kk=0
           ind=1+ii+2*jj+4*kk
           iskip=ncoarse+(ind-1)*ngridmax
           ind_cell=iskip+igrid
           igrid_old = igrid
           igrid=son(ind_cell)
           if(igrid==0.or.j==ilevel)  exit
        end do

        do idim=1,ndim
           xtcell(i,idim)=xg(igrid_old,idim)+xc(ind,idim)
           xtcell(i,idim)=(xtcell(i,idim)-skip_loc(idim))
        end do

        cell_index(i)=ind_cell
        cell_levl(i)=j
     endif
  end do

  
end subroutine get_cell_index4


! CC 03/17
! Get the typical height scale z0 of the gas (mass_all does not work with amr)
subroutine get_height_scale(ilevel,z0,mass_all)
   use amr_commons
   use hydro_commons
   implicit none
   include 'mpif.h'
   integer::ilevel
   integer::i,ind,ncache,igrid,iskip
   integer::info,nleaf,ngrid,nx_loc
   real(dp),dimension(1:3)::skip_loc
   integer::rang,nb_procs,rang_z
   integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_leaf
   real(dp)::dx,vol,scale,z0
   real(kind=8)::mass_loc,mass_max,mass_diff,mass_all
   real(kind=8),allocatable::mass_tab(:)

  nx_loc=icoarse_max-icoarse_min+1
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale
  vol=dx**ndim

  mass_loc=0.0d0

  call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,info)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rang,info)

  allocate(mass_tab(0:nb_procs-1))

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     
     ! Loop over cells
     do ind=1,twotondim        
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=ind_grid(i)+iskip
        end do
        
        ! Gather leaf cells
        nleaf=0
        do i=1,ngrid
           if(son(ind_cell(i))==0)then
              nleaf=nleaf+1
              ind_leaf(nleaf)=ind_cell(i)
           end if
        end do

        do i=1,nleaf
           mass_loc=mass_loc+uold(ind_leaf(i),1)*vol
        end do
     
      end do
 
   end do

   !mass_all=0d0
   !z0=0d0
   !mass_tab(:)=0d0
   !mass_max=0d0
   !nb_procs=0
   !rang=0

   call MPI_ALLGATHER(mass_loc,1,MPI_DOUBLE_PRECISION,mass_tab,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,info)
   call MPI_ALLREDUCE(mass_loc,mass_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
   call MPI_ALLREDUCE(mass_loc,mass_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
        
        ! The height z0 is such that M(z0)=Mmax/2**1.5
   mass_tab(:)=abs(mass_tab(:)-mass_max/(2**1.5)) 
 
   mass_diff=mass_tab(0)
   rang_z=0
   do i=1,nb_procs-1
      if(mass_tab(i)<mass_diff) then  
        mass_diff=mass_tab(i)
        rang_z=i
      end if
   end do

   if(rang==rang_z) then
      z0=abs((xg(ind_grid(1),3)-skip_loc(3)-0.5))*scale
   end if
   call MPI_BCAST(z0,1,MPI_DOUBLE_PRECISION,rang_z,MPI_COMM_WORLD,info)

end subroutine get_height_scale
        
