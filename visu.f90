!################################################################################
!This file is part of Incompact3d.
!
!Incompact3d
!Copyright (c) 2012 Eric Lamballais and Sylvain Laizet
!eric.lamballais@univ-poitiers.fr / sylvain.laizet@gmail.com
!
!    Incompact3d is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation.
!
!    Incompact3d is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with the code.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!    We kindly request that you cite Incompact3d in your publications and 
!    presentations. The following citations are suggested:
!
!    1-Laizet S. & Lamballais E., 2009, High-order compact schemes for 
!    incompressible flows: a simple and efficient method with the quasi-spectral 
!    accuracy, J. Comp. Phys.,  vol 228 (15), pp 5989-6015
!
!    2-Laizet S. & Li N., 2011, Incompact3d: a powerful tool to tackle turbulence 
!    problems with up to 0(10^5) computational cores, Int. J. of Numerical 
!    Methods in Fluids, vol 67 (11), pp 1735-1757
!################################################################################

!############################################################################
!
subroutine VISU_INSTA (ux1,uy1,uz1,phi1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
     ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
     ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG,uvisu)   
!
!############################################################################

USE param
USE variables
USE decomp_2d
USE decomp_2d_io

implicit none

TYPE(DECOMP_INFO) :: phG
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,phi1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3
real(mytype),dimension(xszV(1),xszV(2),xszV(3)) :: uvisu 

integer :: code,icomplet
integer :: ijk,nvect1,nvect2,nvect3,i,j,k
character(len=20) nfichier,nfichier1
character(len=20) :: filename

!FZ_P=10. !!! **************omar

nvect1=xsize(1)*xsize(2)*xsize(3)
!x-derivatives
call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
!y-derivatives
call transpose_x_to_y(ux1,td2)
call transpose_x_to_y(uy1,te2)
call transpose_x_to_y(uz1,tf2)
call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
!!z-derivatives
call transpose_y_to_z(td2,td3)
call transpose_y_to_z(te2,te3)
call transpose_y_to_z(tf2,tf3)
call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
!!all back to x-pencils
call transpose_z_to_y(ta3,td2)
call transpose_z_to_y(tb3,te2)
call transpose_z_to_y(tc3,tf2)
call transpose_y_to_x(td2,tg1)
call transpose_y_to_x(te2,th1)
call transpose_y_to_x(tf2,ti1)
call transpose_y_to_x(ta2,td1)
call transpose_y_to_x(tb2,te1)
call transpose_y_to_x(tc2,tf1)
!du/dx=ta1 du/dy=td1 and du/dz=tg1
!dv/dx=tb1 dv/dy=te1 and dv/dz=th1
!dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1


!############################################################################
!VORTICITY
di1=0.
do ijk=1,nvect1
   di1(ijk,1,1)=sqrt((tf1(ijk,1,1)-th1(ijk,1,1))**2+&
        (tg1(ijk,1,1)-tc1(ijk,1,1))**2+&
        (tb1(ijk,1,1)-td1(ijk,1,1))**2)
enddo
uvisu=0.
call fine_to_coarseV(1,di1,uvisu)
990 format('vort',I3.3)
write(filename, 990) 1
call decomp_2d_write_one(1,uvisu,filename,2)
!call decomp_2d_write_one(nx_global,ny_global,nz_global,&
!     1,di1,filename)
!############################################################################

!############################################################################
!VELOCITY
uvisu=0.
call fine_to_coarseV(1,ux1,uvisu)
993 format('ux',I3.3)
      write(filename, 993) 1
call decomp_2d_write_one(1,uvisu,filename,2)
!call decomp_2d_write_one(nx_global,ny_global,nz_global,&
!           1,ux1,filename)
uvisu=0.
call fine_to_coarseV(1,uy1,uvisu)
994 format('uy',I3.3)
      write(filename, 994) 1
call decomp_2d_write_one(1,uvisu,filename,2)
!call decomp_2d_write_one(nx_global,ny_global,nz_global,&
!           1,uy1,filename)
uvisu=0.
call fine_to_coarseV(1,uz1,uvisu)
995 format('uz',I3.3)
      write(filename, 995) 1
call decomp_2d_write_one(1,uvisu,filename,2)
!call decomp_2d_write_one(nx_global,ny_global,nz_global,&
!           1,uz1,filename)
  
!############################################################################

!############################################################################
!PASSIVE SCALAR
if (iscalar==1) then
uvisu=0.
call fine_to_coarseV(1,phi1,uvisu)
996 format('phi',I3.3)
   write(filename, 996) itime/imodulo
   call decomp_2d_write_one(1,uvisu,filename,2)
!   call decomp_2d_write_one(nx_global,ny_global,nz_global,&
!        1,phi1,filename)
endif
!############################################################################

!############################################################################
!PRESSURE
!IT IS IN A SEPARATE SUBROUTINE
!############################################################################
end subroutine VISU_INSTA

!############################################################################
!
subroutine STATISTIC(ux1,uy1,uz1,phi1,ta1,umean,umean_temp,vmean,wmean,phimean,uumean,vvmean,wwmean,&
     uvmean,uwmean,vwmean,phiphimean,tmean)
!
!############################################################################

USE param
USE variables
USE decomp_2d
USE decomp_2d_io

implicit none

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,phi1
real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: umean,umean_temp,vmean,wmean,uumean,vvmean,wwmean,uvmean,uwmean,vwmean,tmean
real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: phimean, phiphimean
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1

! mean pressure

!umean=ux1
!if (mod(itime,imodulo)==0) then
!umean=0.
!endif
!vmean=0.
!wmean=0.
!phimean=0. 
!uumean=0.
!vvmean=0.
!wwmean=0.
!uvmean=0. 
call fine_to_coarseS(1,ux1,tmean)
umean(:,:,:)=umean(:,:,:)+tmean(:,:,:)

call fine_to_coarseS(1,ux1,tmean) ! omar
umean_temp(:,:,:)=umean_temp(:,:,:)+tmean(:,:,:) 

!vmean=uy1
call fine_to_coarseS(1,uy1,tmean)
vmean(:,:,:)=vmean(:,:,:)+tmean(:,:,:)

!wmean=uz1
call fine_to_coarseS(1,uz1,tmean)
wmean(:,:,:)=wmean(:,:,:)+tmean(:,:,:)

if (iscalar==1) then
   !phimean=phi1
   call fine_to_coarseS(1,phi1,tmean)
   phimean(:,:,:)=phimean(:,:,:)+tmean(:,:,:)
endif

!uumean=ux1*ux1
ta1(:,:,:)=ux1(:,:,:)*ux1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
uumean(:,:,:)=uumean(:,:,:)+tmean(:,:,:)

!vvmean=uy1*uy1
ta1(:,:,:)=uy1(:,:,:)*uy1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
vvmean(:,:,:)=vvmean(:,:,:)+tmean(:,:,:)

!wwmean=uz1*uz1
ta1(:,:,:)=uz1(:,:,:)*uz1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
wwmean(:,:,:)=wwmean(:,:,:)+tmean(:,:,:)

!uvmean=ux1*uy1
ta1(:,:,:)=ux1(:,:,:)*uy1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
uvmean(:,:,:)=uvmean(:,:,:)+tmean(:,:,:)

!uwmean=ux1*uz1
ta1(:,:,:)=ux1(:,:,:)*uz1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
uwmean(:,:,:)=uwmean(:,:,:)+tmean(:,:,:)

!vwmean=uy1*uz1
ta1(:,:,:)=uy1(:,:,:)*uz1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
vwmean(:,:,:)=vwmean(:,:,:)+tmean(:,:,:)

if (iscalar==1) then
   !phiphimean=phi1*phi1
   ta1(:,:,:)=phi1(:,:,:)*phi1(:,:,:)
   call fine_to_coarseS(1,ta1,tmean)
   phiphimean(:,:,:)=phiphimean(:,:,:)+tmean(:,:,:)
endif

!for a verification
!call decomp_2d_write_one(nx_global,ny_global,nz_global,&
!           1,ta1,'compa.dat')

if (mod(itime,imodulo)==0) then

   call decomp_2d_write_one(1,umean,'umean.dat',1)
   call decomp_2d_write_one(1,umean_temp,'umean_temp.dat',1)  ! omar
   call decomp_2d_write_one(1,vmean,'vmean.dat',1)
   call decomp_2d_write_one(1,wmean,'wmean.dat',1)
   call decomp_2d_write_one(1,uumean,'uumean.dat',1)
   call decomp_2d_write_one(1,vvmean,'vvmean.dat',1)
   call decomp_2d_write_one(1,wwmean,'wwmean.dat',1)
   call decomp_2d_write_one(1,uvmean,'uvmean.dat',1)
   call decomp_2d_write_one(1,uwmean,'uwmean.dat',1)
   call decomp_2d_write_one(1,vwmean,'vwmean.dat',1)
   if (nrank==0) print *,'write stat arrays velocity done!'
   if (iscalar==1) then
      call decomp_2d_write_one(1,phimean,'phimean.dat',1)
      call decomp_2d_write_one(1,phiphimean,'phiphimean.dat',1)
      if (nrank==0) print *,'write stat arrays scalar done!'
   endif
!   call decomp_2d_write_one(nx_global,ny_global,nz_global,1,ux1,'compa.dat')
umean_temp=0. ! omar  
!umean=0.
!wmean=0.
!uumean=0.
!vvmean=0.
!wwmean=0.
!uvmean=0.
!uwmean=0.
!vwmean=0.
endif

end subroutine STATISTIC


!############################################################################ omar

subroutine vort_stats (ux1,uy1,uz1,td1,te1,tf1,tg1,th1,ti1,&
     ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3,vort_mean,vort_2_mean)   !!! OMAR
!
!############################################################################ omar

USE param
USE variables
USE decomp_2d
USE decomp_2d_io
implicit none

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,vort,vort_2
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,di2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,di3
real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: vort_mean,vort_2_mean
integer :: ijk,nvect1,nvect2,nvect3,i,j,k
integer :: code,icomplet

character(len=20) nfichier,nfichier1
character(len=20) :: filename

nvect1=xsize(1)*xsize(2)*xsize(3)
!x-derivatives
call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
!y-derivatives
call transpose_x_to_y(ux1,td2)!
call transpose_x_to_y(uy1,te2)!
call transpose_x_to_y(uz1,tf2)!
call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)!
call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)!
call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)!
!!z-derivatives
call transpose_y_to_z(td2,td3)!
call transpose_y_to_z(te2,te3)!
call transpose_y_to_z(tf2,tf3)
call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
!!all back to x-pencils
call transpose_z_to_y(ta3,td2)!
call transpose_z_to_y(tb3,te2)
call transpose_z_to_y(tc3,tf2)
call transpose_y_to_x(td2,tg1)!
call transpose_y_to_x(te2,th1)!
call transpose_y_to_x(tf2,ti1)
call transpose_y_to_x(ta2,td1)
call transpose_y_to_x(tb2,te1)!
call transpose_y_to_x(tc2,tf1)!

!VORTICITY
di1=0.
vort_2=0.
do ijk=1,nvect1
   di1(ijk,1,1)=sqrt((tf1(ijk,1,1)-th1(ijk,1,1))**2+&
        (tg1(ijk,1,1)-tc1(ijk,1,1))**2+&
        (tb1(ijk,1,1)-td1(ijk,1,1))**2)
enddo

vort_mean(:,:,:)=vort_mean(:,:,:)+di1(:,:,:)


vort_2(:,:,:)=di1(:,:,:)*di1(:,:,:)
vort_2_mean(:,:,:)=vort_2_mean(:,:,:)+vort_2(:,:,:)

if (mod(itime,imodulo)==0) then
   call decomp_2d_write_one(1,vort_mean,'vort_mean.dat',1)
   call decomp_2d_write_one(1,vort_2_mean,'vort_2_mean.dat',1)
   if (nrank==0) print *,'write stat arrays VORTICITY done!'     
endif

end subroutine  vort_stats ! omar

!############################################################################ omar

subroutine pp_stats (pp3,ta1,tb1,di1,ta2,tb2,di2,&
     ta3,di3,nxmsize,nymsize,nzmsize,phG,ph2,ph3,uvisu,p_mean,p_2_mean)   !!! OMAR
!
!############################################################################ omar


USE param
USE variables
USE decomp_2d
USE decomp_2d_io

implicit none

integer :: nxmsize,nymsize,nzmsize
TYPE(DECOMP_INFO) :: phG,ph2,ph3
real(mytype),dimension(xszV(1),xszV(2),xszV(3)) :: uvisu 

real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),nzmsize) :: pp3 
!Z PENCILS NXM NYM NZM-->NXM NYM NZ
real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),zsize(3)) :: ta3,di3
!Y PENCILS NXM NYM NZ -->NXM NY NZ
real(mytype),dimension(ph3%yst(1):ph3%yen(1),nymsize,ysize(3)) :: ta2
real(mytype),dimension(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)) :: tb2,di2
!X PENCILS NXM NY NZ  -->NX NY NZ
real(mytype),dimension(nxmsize,xsize(2),xsize(3)) :: ta1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tb1,di1 ,p_2

integer :: code,icomplet
integer :: ijk,nvect1,nvect2,nvect3,i,j,k
real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: p_mean,p_2_mean

!WORK Z-PENCILS
call interiz6(ta3,pp3,di3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
     (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
!WORK Y-PENCILS
call transpose_z_to_y(ta3,ta2,ph3) !nxm nym nz
call interiy6(tb2,ta2,di2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
     (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
!WORK X-PENCILS
call transpose_y_to_x(tb2,ta1,ph2) !nxm ny nz
call interi6(tb1,ta1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
     nxmsize,xsize(1),xsize(2),xsize(3),1)
!The pressure field on the main mesh is in tb1
!PRESSURE
uvisu=0.
call fine_to_coarseV(1,tb1,uvisu)
p_mean(:,:,:)=p_mean(:,:,:)+tb1(:,:,:)

!p_2(:,:,:)=tb1(:,:,:)*tb1(:,:,:)
!uvisu=0.
!call fine_to_coarseV(1,p_2,uvisu)
p_2_mean(:,:,:)=p_2_mean(:,:,:)+tb1(:,:,:)*tb1(:,:,:)

if (mod(itime,imodulo)==0) then
   call decomp_2d_write_one(1,p_mean,'p_mean',2)
   call decomp_2d_write_one(1,p_2_mean,'p_2_mean',2)
   if (nrank==0) print *,'write stat arrays PRESURE done!'     
endif

end subroutine pp_stats


!############################################################################
!
subroutine VISU_PRE (pp3,ta1,tb1,di1,ta2,tb2,di2,&
     ta3,di3,nxmsize,nymsize,nzmsize,phG,ph2,ph3,uvisu)
!
!############################################################################

USE param
USE variables
USE decomp_2d
USE decomp_2d_io

implicit none

integer :: nxmsize,nymsize,nzmsize
TYPE(DECOMP_INFO) :: phG,ph2,ph3
real(mytype),dimension(xszV(1),xszV(2),xszV(3)) :: uvisu 

real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),nzmsize) :: pp3 
!Z PENCILS NXM NYM NZM-->NXM NYM NZ
real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),zsize(3)) :: ta3,di3
!Y PENCILS NXM NYM NZ -->NXM NY NZ
real(mytype),dimension(ph3%yst(1):ph3%yen(1),nymsize,ysize(3)) :: ta2
real(mytype),dimension(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)) :: tb2,di2
!X PENCILS NXM NY NZ  -->NX NY NZ
real(mytype),dimension(nxmsize,xsize(2),xsize(3)) :: ta1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tb1,di1 

integer :: code,icomplet
integer :: ijk,nvect1,nvect2,nvect3,i,j,k
character(len=20) nfichier,nfichier1
character(len=20) :: filename

!WORK Z-PENCILS
call interiz6(ta3,pp3,di3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
     (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
!WORK Y-PENCILS
call transpose_z_to_y(ta3,ta2,ph3) !nxm nym nz
call interiy6(tb2,ta2,di2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
     (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
!WORK X-PENCILS
call transpose_y_to_x(tb2,ta1,ph2) !nxm ny nz
call interi6(tb1,ta1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
     nxmsize,xsize(1),xsize(2),xsize(3),1)
!The pressure field on the main mesh is in tb1
!PRESSURE
uvisu=0.
call fine_to_coarseV(1,tb1,uvisu)
990 format('pp',I3.3)
      write(filename, 990) 1
call decomp_2d_write_one(1,uvisu,filename,2)
!call decomp_2d_write_one(nx_global,ny_global,nz_global,&
!           1,tb1,filename)

end subroutine VISU_PRE

!############################################################################
!
subroutine visu_paraview()
!
!############################################################################
USE param

USE variables
  implicit none
  integer(4) :: nfiles, icrfile, file1, filen, ifile, dig1, dig2, dig3, dig4
  real(4), allocatable :: y1(:),y3(:)
  integer(4) :: i, j, k, num, aig, ii, nfil

!IF THE DATA ARE STORED WITH 3 DIGITS, IE UX001,UX002,ETC.
  character(3) :: chits
!IF THE DATA ARE STORED WITH 4 DIGITS, IE UX0001,UX0002,ETC.
!  character(4) :: chits

  !write (*,*) 'nx, ny, nz   - Incompact3D'
  !read (*,*) nx, ny, nz
  !write (*,*) 'xlx, yly, zlz   - Incompact3D'
  !read (*,*) xlx, yly, zlz
  !write (*,*) 'nclx, ncly, nclz   - Incompact3D'
 ! read (*,*) nclx, ncly, nclz
  !write (*,*) 'n files, first file, last file'
  !read (*,*) nfiles,file1, filen
  !write (*,*) 'Stretching in the y direction (Y=1/N=0)?'
  !read (*,*) istret
 
  
  nfiles=ilast/imodulo
  file1=1
  filen=nfiles

  istret=1

  if (nclx==0) dx=xlx/nx
 if (nclx==1 .or. nclx==2) dx=xlx/(nx-1.)
  if (ncly==0) dy=yly/ny
  if (ncly==1.or.ncly==2) dy=yly/(ny-1.)
 if (nclz==0) dz=zlz/nz
  if (nclz==1.or.nclz==2) dz=zlz/(nz-1.)
  dt=1.

  allocate(y1(nx))
  allocate(y3(nz))
  do i=1,nx
     y1(i)=(i-1)*dx
  enddo
  if (istret==1) then
     print *,'We need to read the yp.dat file'
     open(12,file='yp.dat',form='formatted',status='unknown')
     do j=1,ny
        read(12,*) yp(j)
     enddo
     close(12)
  else
     do j=1,ny
        yp(j)=(j-1)*dy
     enddo
  endif
  do k=1,nz
     y3(k)=(k-1)*dz
  enddo


  nfil=41
  open(nfil,file='visu.xdmf')

  write(nfil,'(A22)')'<?xml version="1.0" ?>'
  write(nfil,*)'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
  write(nfil,*)'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
  write(nfil,*)'<Domain>'
  write(nfil,*)'    <Topology name="topo" TopologyType="3DRectMesh"'
  write(nfil,*)'        Dimensions="',nz,ny,nx,'">'
  write(nfil,*)'    </Topology>'
  write(nfil,*)'    <Geometry name="geo" Type="VXVYVZ">'
  write(nfil,*)'    <DataItem Dimensions="',nx,'" NumberType="Float" Precision="4" Format="XML">'
  write(nfil,*)'    ',y1(:) 
  write(nfil,*)'    </DataItem>'
  write(nfil,*)'    <DataItem Dimensions="',ny,'" NumberType="Float" Precision="4" Format="XML">'
  write(nfil,*)'    ',yp(:) 
  write(nfil,*)'    </DataItem>'
  write(nfil,*)'    <DataItem Dimensions="',nz,'" NumberType="Float" Precision="4" Format="XML">'
  write(nfil,*)'    ',y3(:) 
  write(nfil,*)'    </DataItem>'
  write(nfil,*)'    </Geometry>'
  write(nfil,'(/)')
  write(nfil,*)'    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
  write(nfil,*)'        <Time TimeType="HyperSlab">'
  write(nfil,*)'            <DataItem Format="XML" NumberType="Float" Dimensions="3">'
  write(nfil,*)'           <!--Start, Stride, Count-->'
  write(nfil,*)'            0.0',dt
  write(nfil,*)'            </DataItem>'
  write(nfil,*)'        </Time>'

  do ifile = file1, filen

!IF THE DATA ARE STORED WITH 4 DIGITS, IE UX0001,UX0002,ETC.
  !   dig1 =   ifile/1000 + 48
  !   dig2 = ( ifile - 1000*( ifile/1000 ) )/100 + 48
  !   dig3 = ( ifile - 100*( ifile/100 ) )/10 + 48
  !   dig4 = ( ifile - 10*( ifile/10 ) )/1 + 48
  !   chits(1:4) = char(dig1)//char(dig2)//char(dig3)//char(dig4)    

!IF THE DATA ARE STORED WITH 3 DIGITS, IE UX001,UX002,ETC.
    dig1 =   ifile/100 + 48
    dig2 = ( ifile - 100*( ifile/100 ) )/10 + 48
    dig3 = ( ifile - 10*( ifile/10 ) )/1 + 48
    chits(1:3) = char(dig1)//char(dig2)//char(dig3)

     write(*,*) ifile, 'file'//chits

     write(nfil,'(/)')
     write(nfil,*)'        <Grid Name="'//chits//'" GridType="Uniform">'
     write(nfil,*)'            <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
     write(nfil,*)'            <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'
!SINGLE PRECISION-->Precision=4
!DOUBLE PRECISION-->Precision=8
     write(nfil,*)'            <Attribute Name="ux" Center="Node">'
     write(nfil,*)'               <DataItem Format="Binary" '
     write(nfil,*)'                DataType="Float" Precision="8" Endian="little"'
     write(nfil,*)'                Dimensions="',nz,ny,nx,'">'
     write(nfil,*)'                  ux'//chits
     write(nfil,*)'               </DataItem>'
     write(nfil,*)'            </Attribute>'

!it is possible to add as much field as you want for example uy

    write(nfil,*)'            <Attribute Name="uy" Center="Node">'
    write(nfil,*)'               <DataItem Format="Binary" '
    write(nfil,*)'                DataType="Float" Precision="8" Endian="little"'
    write(nfil,*)'                Dimensions="',nz,ny,nx,'">'
    write(nfil,*)'                  uy'//chits
    write(nfil,*)'               </DataItem>'
    write(nfil,*)'            </Attribute>'

    write(nfil,*)'            <Attribute Name="uz" Center="Node">'
    write(nfil,*)'               <DataItem Format="Binary" '
    write(nfil,*)'                DataType="Float" Precision="8" Endian="little"'
    write(nfil,*)'                Dimensions="',nz,ny,nx,'">'
    write(nfil,*)'                  uz'//chits
    write(nfil,*)'               </DataItem>'
    write(nfil,*)'            </Attribute>'


    write(nfil,*)'            <Attribute Name="vort" Center="Node">'
    write(nfil,*)'               <DataItem Format="Binary" '
    write(nfil,*)'                DataType="Float" Precision="8" Endian="little"'
    write(nfil,*)'                Dimensions="',nz,ny,nx,'">'
    write(nfil,*)'                  vort'//chits
    write(nfil,*)'               </DataItem>'
    write(nfil,*)'            </Attribute>'

    write(nfil,*)'            <Attribute Name="pp" Center="Node">'
    write(nfil,*)'               <DataItem Format="Binary" '
    write(nfil,*)'                DataType="Float" Precision="8" Endian="little"'
    write(nfil,*)'                Dimensions="',nz,ny,nx,'">'
    write(nfil,*)'                  pp'//chits
    write(nfil,*)'               </DataItem>'
    write(nfil,*)'            </Attribute>'

     write(nfil,*)'        </Grid>'

  enddo
  write(nfil,'(/)')
  write(nfil,*)'    </Grid>'
  write(nfil,*)'</Domain>'
  write(nfil,'(A7)')'</Xdmf>'
  close(nfil)

end subroutine visu_paraview


!############################################################################
!
subroutine verif(Re,itime,ifirst,nx1, ny1, nz1)
!
!############################################################################
  implicit none

  integer :: nx1, ny1, nz1,xitime,itime,ifirst
  real(8),dimension(nx1,ny1,nz1) :: umean,uumean
  integer :: i,j,k,count,nfil
  real(8) :: pi,x,y,u_to,Re,Re_t_up,Re_t_low
  real(8),dimension(ny1) :: yp,ypi,qstat,qstat2
  real(8),dimension(nx1) :: y1
  real(8),dimension(ny1) :: y2
  real(8),dimension(nz1) :: y3
  logical :: exist

  character(10) :: time  
  integer,dimension(8) :: values  
  call date_and_time(TIME=time)
  call date_and_time(VALUES=values)
  print*,'== Machine Time == '
  print '(8i5)', values

pi=acos(-1.)
u_to=180./Re 

 xitime=itime-ifirst+1

   open (15,file='yp.dat',form='formatted',status='unknown')
   do j=1,ny1
      read(15,*) yp(j)
   enddo
   close(15)


!! U Perturbation statistics !!

 OPEN(11,FILE='umean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8)
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) umean(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
  ENDDO

  CLOSE(11)
OPEN(11,FILE='uumean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8)
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) uumean(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
  ENDDO
  CLOSE(11)

umean=umean/xitime
uumean=uumean/xitime


uumean=(uumean-umean*umean)!/u_to/u_to
qstat=0.
qstat2=0.
         do k=1,nz1
         do j=1,ny1
         do i=1,nx1
            qstat(j)=qstat(j)+uumean(i,j,k)
            qstat2(j)=qstat2(j)+umean(i,j,k)
         enddo
         enddo
         enddo

         qstat(:)=qstat(:)/nx1/nz1
	 qstat2(:)=qstat2(:)/nx1/nz1

  Re_t_up=sqrt((qstat2(ny1-1)-qstat2(ny1))/(yp(ny1)-yp(ny1-1))/Re)*Re
  Re_t_low=sqrt((qstat2(2)-qstat2(1))/(yp(2)-yp(1))/Re)*Re
  print *,itime,xitime,Re_t_up,Re_t_low,abs(Re_t_up-Re_t_low)
        
u_to=sqrt((qstat2(ny1-1)-qstat2(ny1))/(yp(ny1)-yp(ny1-1))/Re)*Re+sqrt((qstat2(2)-qstat2(1))/(yp(2)-yp(1))/Re)*Re
u_to=u_to*0.5/Re

print *,'U_TO',u_to
  open(10,file='upert.dat',status='unknown',form='formatted')
  do j=1,ny1
     write(10,*) yp(j)*u_to*Re,qstat(j)/u_to/u_to,qstat2/u_to/u_to
  enddo
  close(10)

!! V Pertubation Statistics !!

 OPEN(11,FILE='vmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8)
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) umean(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
  ENDDO
  CLOSE(11)
OPEN(11,FILE='vvmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8)
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) uumean(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
  ENDDO
  CLOSE(11)

umean=umean/xitime
uumean=uumean/xitime

uumean=(uumean-umean*umean)!/u_to/u_to

         do k=1,nz1
         do j=1,ny1
         do i=1,nx1
            qstat(j)=qstat(j)+uumean(i,j,k)
         enddo
         enddo
         enddo

         qstat(:)=qstat(:)/nx1/nz1


  open(10,file='vpert.dat',status='unknown',form='formatted')
  do j=1,ny1
     write(10,*) yp(j)*u_to*Re,qstat(j)/u_to/u_to
  enddo
  close(10)

!! W Pertubation Statistics !!

 OPEN(11,FILE='wmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8)
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) umean(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
  ENDDO
  CLOSE(11)
OPEN(11,FILE='wwmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8)
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) uumean(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
  ENDDO
  CLOSE(11)

umean=umean/xitime
uumean=uumean/xitime

uumean=(uumean-umean*umean)!/u_to/u_to

         do k=1,nz1
         do j=1,ny1
         do i=1,nx1
            qstat(j)=qstat(j)+uumean(i,j,k)
         enddo
         enddo
         enddo

         qstat(:)=qstat(:)/nx1/nz1

  open(10,file='wpert.dat',status='unknown',form='formatted')
  do j=1,ny1
     write(10,*) yp(j)*u_to*Re,qstat(j)/u_to/u_to
  enddo
  close(10)


!************* write shear velocities and Re_t at each call of the subroutine ************
  
  inquire(file="time_stats.dat", exist=exist)
  if (exist) then
    open(12, file="time_stats.dat", status="old", position="append", action="write")
  else
    open(12, file="time_stats.dat", status="new", action="write")
  end if
  write(12, *) itime,xitime,u_to,Re_t_up,Re_t_low,abs(Re_t_up-Re_t_low)
  close(12)
!******************************************************************************************
    end subroutine verif
!
