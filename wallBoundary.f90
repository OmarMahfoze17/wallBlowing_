!############################################################################
!
subroutine wallBoundary (ux,uy,uz)
    !
    !############################################################################

    USE decomp_2d
    USE variables
    USE param
    USE var
    USE MPI


    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    !real(mytype),dimension(nx,nz)::uxWall,uyWall,uzWall
    integer :: i,j,k,ii,kk,code,coun,nx_Act,nz_Act
    real :: uxWallAmpl,uyWallAmpl,uzWallAmpl,lx0,lx1,lz0,lz1,lx01,lz01
    real :: alpha_x,alpha_z,waveLenght_x,waveLenght_z
    integer, dimension(2) :: dims, dummy_coords
    logical, dimension(2) :: dummy_periods

    !****************************************************
    !--------- Define the whole values of the wall ------
    !****************************************************

    !print *, xlx ,zlz,nx,nz,dx,dz
    uyWallAmpl=1.
    nx_Act=2
    nz_Act=5
    alpha_x=0.5
    alpha_z=0.5
    waveLenght_x=xlx/nx_Act
    waveLenght_z=zlz/nz_Act


    lx0=waveLenght_x*(1.-alpha_x)   ! the length of the inactive part of the wall measured in X-direction
    lx1=waveLenght_x*alpha_x	 ! the length of the active part of the wall measured in X-direction
    lz0=waveLenght_z*(1.-alpha_z)    ! the length of the inactive part of the wall measured in z-direction
    lz1=waveLenght_z*alpha_z	 ! the length of the active part of the wall measured in z-direction
    print *, waveLenght_x,waveLenght_z,lx0,lx1
    ! determine the processor grid in use
    call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, &
        dims, dummy_periods, dummy_coords, code)


    if (ncly==2) then
        !find j=1 and j=ny
        if (xstart(2)==1) then  !! the Cores next to the lower wall

            kk=0
            do k=1,nz
            ii=0
            do i=1,nx
            if (dx*ii<=lx1 .and. dz*kk<=lz1 .and. k>=xstart(3) .and. k<=xend(3)) then
                uy(i,1,k-xstart(3)+1)=uyWallAmpl        
            endif
            ii=ii+1
            if (dx*ii>lx1+lx0) then
                ii=0
            endif

            enddo    
            kk=kk+1
            if (dz*kk>lz1+lz0) then
                kk=0
            endif
            enddo
        endif
    

    if (ny-(nym/dims(1))==xstart(2)) then  !! the Cores next to the upper wall
        kk=0
        do k=1,nz
        ii=0
        do i=1,nx
        if (dx*ii<=lx1 .and. dz*kk<=lz1 .and. k>=xstart(3) .and. k<=xend(3)) then
            uy(i,xsize(2),k-xstart(3)+1)=-uyWallAmpl        
        endif
        ii=ii+1
        if (dx*ii>lx1+lx0) then
            ii=0
        endif

        enddo    
        kk=kk+1
        if (dz*kk>lz1+lz0) then
            kk=0
        endif
        enddo
    endif
    endif
    return
end subroutine wallBoundary

