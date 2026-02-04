!
!******************************************************************
!
!  Maxwellian Particle Generator
!  Generates np monospecies particles in random positions with a 
!  random velocity that follows the Maxwellian speed distribution
!  in a 3D space
!
!  Variables preceded by 't' means temporary
!
!  Author: Andrew Hong (02/21/2018)
!
!******************************************************************
!
!
	Program DES_Particle_Generator

	implicit none

 	integer i, j, k, sz
	integer, dimension(:), allocatable :: iseed
	integer mmax
	integer time(8)

	double precision :: ub, vb, wb
	! from input file
	double precision, dimension(:), allocatable :: eps, radius, density
	! precomputed values
	double precision, dimension(:), allocatable :: rad1, vol
	integer, dimension(:), allocatable :: np
	

	real*8, dimension(:,:,:), allocatable :: u
	real*8, dimension(:,:), allocatable :: x,y,z
	real*8, dimension(3) :: bvel
	real*8 r1, r2
	real*8 u1, u2, u3
	real*8 xp, yp, zp
	double precision :: xl, yl, zl
	real*8 dia, dist, 
	real*8 density
	real*8 mass, vel, p1, c1, grantemp
	real*8 velc, pc
	real*8 pi
	real*8 :: cut = 1e-8

	logical vlog

	open(unit=10, file="Pgen.in", status='old')
	open(unit=20, file="particle_input.dat", status='unknown')

	read (10,*) mmax
	allocate(eps(mmax));allocate(radius(mmax));allocate(density(mmax));allocate(np(mmax))
	allocate(mass(mmax));allocate(vol(mmax))
	read (10,*) (eps(i), i=1,mmax)
	read (10,*) (radius(i), i=1,mmax)
	read (10,*) (density(i), i=1,mmax)
	read (10,*) xl
	read (10,*) yl
 	read (10,*) zl
	read (10,*) (grantemp(i), i=1,mmax)

	pi = 4.0d0*ATAN(1.0d0)

	if(sum(eps) .gt. 1.0) then
		write(*,*) 'Error: Total solids fraction too high'
		stop
	end if

	do lm = 1,mmax
		vol(lm) = 4.0/3.0*pi*radius(lm)**3
		np(lm) = nint((xl*yl&zl*eps(:))/vol(:))
	end do

	tnp = sum(np(:))
	allocate( u(3,tnp) )
	allocate( x(tnp) ); allocate( y(tnp) ); allocate( z(np) )

!       Specifying the initial distribution of the particles

	call random_seed( size = sz )
	allocate(iseed(sz))
	call date_and_time(values=time)
	iseed = (time(8)*2 + time(7)*3 + time(6)*5)*(/(I, I=1, sz)/)
	call random_seed(PUT=iseed)
	do lm = 1,mmax
	pp = np(lm)
	do i = 1, pp
10		continue
		call random_particle(radius,xp,yp,zp,xl,yl,zl)
		x(i) = xp
		y(i) = yp
		z(i) = zp
		dist = 0d0
		do j = 1, i
			if(j.ne.i) then
				dist = sqrt( (x(i)-x(j))**2 + (y(i)-y(j))**2 +&
		    	(z(i)-z(j))**2 )
				if(dist.le.dia) go to 10
			end if    
		end do
	end do
	end do

	mass(:) = vol(:)*density(:)
	mass(:) = 2*grantemp(:)/mass(:)
	vel = sqrt(mass(:))
	c1 = 4.0d0*pi*(1/(pi*mass))**(1.5)
	p1 = c1*vel**2*exp(-vel**2/mass)

	DO i = 1,np	
		vlog = .True.
		DO WHILE(vlog)
		CALL RANDOM_NUMBER(r1)
		CALL RANDOM_NUMBER(r2)
		velc = 5.0d0*vel*r1
		pc = c1*velc**2*exp(-velc**2/mass)
		IF(pc .GT. p1*r2) vlog = .False.
		ENDDO

	call RANDOM_NUMBER(r1)
	call RANDOM_NUMBER(r2)
	r1 = 2*r1-1
	r2 = r2*2*pi
	u1 = sqrt(1-r1**2)*cos(r2)
	u2 = sqrt(1-r1**2)*sin(r2)
	u3 = r1

	u(1,i) = velc*u1
	u(2,i) = velc*u2
	u(3,i) = velc*u3
	end do

	call bulk_vel(ub, np, u(1,:))
	call bulk_vel(vb, np, u(2,:))
	call bulk_vel(wb, np, u(3,:))
	bvel(1) = ub; bvel(2) = vb; bvel(3) = wb
	write(*,*) bvel

	! attempts to shift velocities until bulk momentum is zero
	do i = 1,3
		do while(abs(bvel(i)).gt.cut)
	   	u(i,:) = u(i,:) - ub
			call bulk_vel(ub, np, u(i,:))
			bvel(i) = ub;
		end do
	end do
	write(*,*) bvel
	

	do i = 1,np
	   write(20,12) x(i), y(i), z(i), radius, density, u(1,i), u(2,i), u(3,i)
	end do

 11     FORMAT (6(d10.4,2x))
 12     FORMAT (8(d10.4,2x))

        stop
        end
!---
! Random particle
!---

        subroutine random_particle(rad,xp1,yp1,zp1,xl1,yl1,zl1)

        integer i, ic
        real*8 rad, xp1, yp1, zp1, xl1, yl1, zl1
        real*8 rad1
        real*8 pxy(3)

	   ic = 100000
           do i = 1, ic
             call random_number(pxy)
             xp1 = dble(pxy(1))*xl1
             yp1 = dble(pxy(2))*yl1
             zp1 = dble(pxy(3))*zl1
             rad1 = 1.01*rad
             if((xp1.ge.rad1).and.(xp1.le.xl1-rad1).and.(yp1.ge.rad1).and. &
		(yp1.le.yl1-rad1).and.(zp1.ge.rad1).and.(zp1.le.zl1-rad1)) exit
           end do
           if(i.gt.ic) then
             print *,'not able to place particle'
	     stop
	   end if

        return
        end subroutine random_particle

!---
! Subroutine: BULK_VEL                                                 
!---
      SUBROUTINE bulk_vel(tvel,tnp,tu)
	
      IMPLICIT NONE

      integer :: tnp, ti

      double precision, intent(out) :: tvel
      real*8, dimension(tnp) :: tu

!---
! Local Variables
!---
	
	tvel = 0.0

	DO ti = 1,tnp		
		tvel = tvel + tu(ti)
	ENDDO
	tvel = tvel/tnp

	return

      END SUBROUTINE BULK_VEL


!
