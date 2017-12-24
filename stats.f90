!############################################################################
!
program verif
!
!############################################################################
  implicit none

  integer, parameter :: nx1=128 , ny1=129, nz1=128
  integer xitime,xitime_temp,itime,ifirst,nx, ny, nz
  real(8),dimension(nx1, ny1, nz1) :: umean,uumean,uvmean,umean_temp
  real(8),dimension(ny1,nz1) :: ux_mean,uux_mean,uuy_mean,uuz_mean,uv_mean_2D
  integer :: i,j,k,count,nfil
  real(8) ::pi,x,y,u_to,Re,Re_t_up,Re_t_low,u_to_temp,Re_t_up_temp,Re_t_low_temp
  real(8),dimension(ny1) :: yp,ypi,qstat,qstat2,qstat2_temp,qstat_u,qstat_v,qstat_w
  real(8),dimension(nx1) :: y1
  real(8),dimension(ny1) :: y2
  real(8),dimension(nz1) :: y3
  logical :: exist

  character(10) :: time  
  integer,dimension(8) :: values 
  open (1,file='comu_stat.dat')
      read(1,*) Re,itime,ifirst,nx, ny, nz,xitime_temp
if (nx.ne.nx1.or.ny.ne.ny1.or.nz.ne.nz1) then 
print*, '###############################################################'
print*,'                dimensions miss match'
print*,nx, ny, nz
print*,nx1, ny1, nz1
print*, '###############################################################'
endif
  close(1)

  call date_and_time(TIME=time)
  call date_and_time(VALUES=values)
  print*,'== Machine Time == '
  print '(8i5)', values

pi=acos(-1.)
u_to=180./Re 

xitime=itime-ifirst+1
!xitime_temp=400
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
close(11)
 OPEN(11,FILE='umean_temp.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8)
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) umean_temp(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
  ENDDO

  CLOSE(11)
OPEN(11,FILE='uumean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8)

OPEN(12,FILE='uvmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8)
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) uumean(I,J,K)
           READ(12,REC=COUNT) uvmean(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
  ENDDO
  CLOSE(12)
  CLOSE(11)

umean=umean/xitime
umean_temp=umean_temp/xitime_temp
uumean=uumean/xitime


uumean=(uumean-umean*umean)!/u_to/u_to
qstat=0.
qstat2=0.
qstat2_temp=0.
         do k=1,nz1
         do j=1,ny1
         do i=1,nx1
            qstat(j)=qstat(j)+uumean(i,j,k)
            qstat2(j)=qstat2(j)+umean(i,j,k)
            qstat2_temp(j)=qstat2_temp(j)+umean_temp(i,j,k)   
            uux_mean(j,k)=uux_mean(j,k)+uumean(i,j,k)
	    ux_mean(j,k)=ux_mean(j,k)+umean(i,j,k)
	    uv_mean_2d(j,k)=uv_mean_2d(j,k)+uvmean(i,j,k)
         enddo
         enddo
         enddo

         qstat(:)=qstat(:)/nx1/nz1
	 qstat2(:)=qstat2(:)/nx1/nz1
         qstat2_temp(:)=qstat2_temp(:)/nx1/nz1

	 ux_mean(:,:)=ux_mean(:,:)/nx1
	 uux_mean(:,:)=uux_mean(:,:)/nx1
	 uv_mean_2d(:,:)=uv_mean_2d(:,:)/nx1
         
         qstat_u=qstat

  Re_t_up=sqrt((qstat2(ny1-1)-qstat2(ny1))/(yp(ny1)-yp(ny1-1))/Re)*Re
  Re_t_low=sqrt((qstat2(2)-qstat2(1))/(yp(2)-yp(1))/Re)*Re

  Re_t_up_temp=sqrt((qstat2_temp(ny1-1)-qstat2_temp(ny1))/(yp(ny1)-yp(ny1-1))/Re)*Re
  Re_t_low_temp=sqrt((qstat2_temp(2)-qstat2_temp(1))/(yp(2)-yp(1))/Re)*Re
  print *,itime,xitime,Re_t_up,Re_t_low,abs(Re_t_up-Re_t_low)
        
u_to=sqrt((qstat2(ny1-1)-qstat2(ny1))/(yp(ny1)-yp(ny1-1))/Re)*Re+sqrt((qstat2(2)-qstat2(1))/(yp(2)-yp(1))/Re)*Re
u_to=u_to*0.5/Re

u_to_temp=sqrt((qstat2_temp(ny1-1)-qstat2_temp(ny1))/(yp(ny1)-yp(ny1-1))/Re)*Re+&
sqrt((qstat2_temp(2)-qstat2_temp(1))/(yp(2)-yp(1))/Re)*Re
u_to_temp=u_to_temp*0.5/Re
print *,'U_TO',u_to
  open(10,file='upert.dat',status='unknown',form='formatted')
  do j=1,ny1
     write(10,*) yp(j)*u_to*Re,qstat(j)/u_to/u_to,qstat2(j)/u_to/u_to
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
	    uuy_mean(j,k)=uuy_mean(j,k)+uumean(i,j,k)
         enddo
         enddo
         enddo

         qstat(:)=qstat(:)/nx1/nz1
	 uuy_mean(:,:)=uuy_mean(:,:)/nx1
         qstat_v=qstat
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
	    uuz_mean(j,k)=uuz_mean(j,k)+uumean(i,j,k)
         enddo
         enddo
         enddo

         qstat(:)=qstat(:)/nx1/nz1
	 uuz_mean(:,:)=uuz_mean(:,:)/nx1
         qstat_w=qstat

  open(10,file='wpert.dat',status='unknown',form='formatted')
  do j=1,ny1
     write(10,*) yp(j)*u_to*Re,qstat(j)/u_to/u_to
  enddo
  close(10)

  open(30,file='u_mean.dat',status='unknown',form='formatted')
  do j=1,ny1
	do k=1,nz1
     	     write(30,*) ux_mean(j,k)/u_to/u_to,uux_mean(j,k)/u_to/u_to,uuy_mean(j,k)/u_to/u_to,&
uuz_mean(j,k)/u_to/u_to,uv_mean_2D(j,k)/u_to/u_to
	enddo
  enddo
  close(30)

  open(10,file='uvwpert.dat',status='unknown',form='formatted')
  do j=1,ny1
     write(10,*) yp(j)*u_to*Re,qstat(j)/u_to/u_to,qstat_u(j)/u_to/u_to,qstat_v(j)/u_to/u_to,qstat_w(j)/u_to/u_to
  enddo
  close(10)
!************* write shear velocities and Re_t at each call of the subroutine ************
  
  inquire(file="time_stats.dat", exist=exist)
  if (exist) then
    open(12, file="time_stats.dat", status="old", position="append", action="write")
  else
    open(12, file="time_stats.dat", status="new", action="write")
  end if
  write(12, *) itime,xitime,u_to,Re_t_up,Re_t_low,abs(Re_t_up-Re_t_low),&
u_to_temp,Re_t_up_temp,Re_t_low_temp,abs(Re_t_up_temp-Re_t_low_temp)
  close(12)
!******************************************************************************************
    end program verif
!
