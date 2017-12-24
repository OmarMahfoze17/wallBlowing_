
PROGRAM mpi_write

implicit none
integer,parameter :: nx=128,ny=129,nz=84
real(8),dimension(ny,NZ) :: FZ2_p,FY2_p,FZ2_n,FY2_n
real(8),dimension(nx,ny,nz) :: FX3_P,FZ3_P,FY3_P,FZ3_N,FY3_N
integer:: i,j,k,COUNT

open(1,file='matlab_data2.vrt')

   do K=1,Nz
        Do j=1,NY
	      READ(1,*) FZ2_p(J,K),FY2_p(J,K),FZ2_n(J,K),FY2_n(J,K)	
	enddo
	
   enddo
CLOSE (1)

print*, maxval(FZ2_p)
DO I=1,NX
   DO J=1,NY
	DO K=1,NZ          
           FX3_P(I,J,K)=0.0
	   FY3_P(I,J,K)=FY2_P(J,K)
	   FZ3_P(I,J,K)=FZ2_P(J,K)
  
	ENDDO
   ENDDO
ENDDO

print*, maxval(FZ3_p)
OPEN(4,FILE='FY_p.dat')
close(4)
open(2,file='FZ_p.dat')
close(2)

OPEN(4,FILE='FY_N.dat')
close(4)
open(2,file='FZ_N.dat')
close(2)

OPEN(4,FILE='FY_p.dat',FORM='UNFORMATTED',ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1     
  DO K=1,nz
     DO J=1,ny
        DO I=1,nx
           WRITE(4,REC=COUNT) FY3_P(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
  ENDDO
 CLOSE(4)

OPEN(4,FILE='FY_N.dat',FORM='UNFORMATTED',ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1     
  DO K=1,nz
     DO J=1,ny
        DO I=1,nx
           WRITE(4,REC=COUNT) FY3_N(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
  ENDDO
 CLOSE(4)

OPEN(5,FILE='FZ_p.dat',FORM='UNFORMATTED',ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1

  DO K=1,nz
     DO J=1,ny
        DO I=1,nx
           WRITE(5,REC=COUNT) FZ3_P(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
  ENDDO
 CLOSE(5)

OPEN(5,FILE='FZ_N.dat',FORM='UNFORMATTED',ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1

  DO K=1,nz
     DO J=1,ny
        DO I=1,nx
           WRITE(5,REC=COUNT) FZ3_N(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
  ENDDO
 CLOSE(5)



OPEN(5,FILE='FZ.dat',FORM='UNFORMATTED',ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1

  DO K=1,nz
     DO J=1,ny
        DO I=1,nx
           read(5,REC=COUNT) FZ3_N(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
  ENDDO
 CLOSE(5)


print*, minval(FZ3_N-FZ3_p)

end PROGRAM  mpi_write



