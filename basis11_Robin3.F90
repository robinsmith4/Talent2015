MODULE CONSTANTS

 Double precision, parameter, public :: pi=3.14159274
 Double precision, parameter, public :: hbar=197.326 !MeV c
 Double precision, parameter, public :: me=937 !MeV / c^2

END MODULE CONSTANTS

PROGRAM BASIS
USE CONSTANTS

 IMPLICIT NONE

        Double precision, external 			     :: DELTA
	Double precision, external 			     :: DELTA2
        Double precision, external			     :: MATRIXELEMENT_V
        Double precision, external			     :: QTRANS
	Integer, dimension(:), allocatable	   	     :: nx, ny, nz, spin
	Integer						     :: nmax,max_orbits,x,y,z,s,i,j,k,nqs,qnums
	Integer, dimension(:,:), allocatable	             :: states
	Double precision, dimension(:), allocatable	     :: E
	Integer, dimension(5)				     :: HOLD
        Double precision				     :: Density
	Double precision				     :: L,Vol,r0
        Double precision				     :: ee,a0,prefac,Kf
	Integer						     :: A,gs
	Double precision, dimension(:,:), allocatable	     :: H_hf
    	Double precision, dimension(:,:), allocatable	     :: C_hf
	Double precision			   	     :: add
        Integer						     :: m,n,gam,delt


 nmax=1
 max_orbits = 1300
 qnums = 5
 A=2
 gs=2
 Kf=1.0
 r0=1.0

 ee=1.4399764
 a0=52917.72
! prefac=ee/(2.0*a0)

 prefac=(hbar**2)/(2.*me)

! prefac=1

 write(6,*)prefac

 Density=1e-17 !fm^-3

 Vol = (A/Density)

 L = Vol**(1./3.)

! L=((6*(pi**2)*A)/(gs*Kf**3))**(1./3.)

 write(6,*)L, pi
 

  allocate ( nx(max_orbits), ny(max_orbits), nz(max_orbits), spin(max_orbits) )
  allocate ( E(max_orbits) )
  allocate ( states(qnums,max_orbits) )
  allocate ( H_hf(max_orbits,max_orbits) )
  allocate ( C_hf(A,max_orbits) )

 nx=0
 ny=0
 nz=0
 spin=0


i=0


 do nqs = 0, nmax*nmax*3

  do x=-nmax,nmax
   do y=-nmax,nmax
    do z=-nmax,nmax
     do s =-1,1,2


	 if (x*x+y*y+z*z /=  nqs)  cycle
	 i = i+1
 
	 nx(i)=x
	 ny(i)=y
	 nz(i)=z
	 spin(i)=s
	 E(i)=nx(i)**2 + ny(i)**2 + nz(i)**2

!         write(6,*)E(i)

	 states(1,i) = nx(i)
	 states(2,i) = ny(i)
	 states(3,i) = nz(i)
	 states(4,i) = spin(i)
	 states(5,i) = E(i)

         E(i)=prefac*((2*pi/L)**2) * (nx(i)**2 + ny(i)**2 + nz(i)**2)


     end do
    end do
   end do
  end do

 end do


 DO j=1,i
  write(6,*)j,states(1,j),states(2,j),states(3,j),states(4,j),states(5,j),E(j)
 END DO

! write(6,*)MATRIXELEMENT_V (1,12,7,10,L,states,max_orbits,qnums)

! First guess for C_hf
 do m=1,A
  do n=1,i
   C_hf(m,n)=0
   if(m.eq.n) C_hf(m,n)=1
  end do
!  write(6,*)C_hf(m,1),C_hf(m,2),C_hf(m,3),C_hf(m,4),C_hf(m,5)
 end do

 do m=1,i
  do n=1,i
   add=0
   H_hf(m,n) = E(m)*DELTA(m,n)
   do j=1,A
    do gam=1,i
     do delt=1,i
!      add=add+(C_hf(j,gam)*C_hf(j,delt)*MATRIXELEMENT_V(m,gam,n,delt,L,states,max_orbits,qnums))
!      if (MATRIXELEMENT_V(m,gam,n,delt,L,states,max_orbits,qnums).gt.0.00001) then
!       write(6,*)m,gam,n,delt,MATRIXELEMENT_V(m,gam,n,delt,L,states,max_orbits,qnums)
!      endif
     end do
    end do
   end do
   H_hf(m,n)=H_hf(m,n)+add
!   write(6,*)'Element added ',m,n,H_hf(m,n)
  end do
 end do

 do j=1,14
  do m=1,14
   write(6,*)MATRIXELEMENT_V(m,j,m,j,L,states,max_orbits,qnums)
  enddo
 enddo


     deallocate ( nx, ny, nz, spin, E )




END PROGRAM  BASIS


!-------------------------------------------------------------------------

 FUNCTION DELTA (n,m) RESULT (r)

 Integer n,m
 Double precision r

 r=0
 if ( n.eq.m ) then 
  r=1
 end if
 END FUNCTION DELTA

!-------------------------------------------------------------------------

 FUNCTION QTRANS (n,m,states,max_orbits,qnums,L) RESULT (q2)

 integer n,m,max_orbits,qnums,states(qnums,max_orbits)
 Double precision L,q2,pi
 integer nx(2),ny(2),nz(2)

 pi=3.1415

 nx(1)=states(1,m)
 ny(1)=states(2,m)
 nz(1)=states(3,m)

 nx(2)=states(1,n)
 ny(2)=states(2,n)
 nz(2)=states(3,n)

 q2=((2*pi/L)**2)*( (nx(2)-nx(1))**2 + (ny(2)-ny(1))**2 + (nz(2)-nz(1))**2 )

 END FUNCTION QTRANS

!-------------------------------------------------------------------------

 FUNCTION MATRIXELEMENT_V (p,q,r,s,L,states,max_orbits,qnums) RESULT (V)
 USE constants
 integer p,q,r,s,max_orbits,qnums
 Double precision L,Qrp,Qsp,DELTA
 Double precision V1,V2,V3,V4,V
 integer states(qnums,max_orbits)

 integer nx(4),ny(4),nz(4),sz(4)

 nx(1)=states(1,p)
 ny(1)=states(2,p)
 nz(1)=states(3,p)
 sz(1)=states(4,p)

 nx(2)=states(1,q)
 ny(2)=states(2,q)
 nz(2)=states(3,q)
 sz(2)=states(4,q)

 nx(3)=states(1,r)
 ny(3)=states(2,r)
 nz(3)=states(3,r)
 sz(3)=states(4,r)

 nx(4)=states(1,s)
 ny(4)=states(2,s)
 nz(4)=states(3,s)
 sz(4)=states(4,s)

 
 Qrp = ((2*pi/L)**2)*( (nx(3)-nx(1))**2 + (ny(3)-ny(1))**2 + (nz(3)-nz(1))**2 )
 Qsp = ((2*pi/L)**2)*( (nx(4)-nx(1))**2 + (ny(4)-ny(1))**2 + (nz(4)-nz(1))**2 )


 V1 = (4*pi/L**3)*DELTA(nx(1)+nx(2),nx(3)+nx(4)) * DELTA(ny(1)+ny(2),ny(3)+ny(4)) * DELTA(nz(1)+nz(2),nz(3)+nz(4))

 V2 = DELTA(sz(1),sz(3)) * DELTA(sz(2),sz(4)) * (1 - DELTA(nx(1),nx(3)) * DELTA(ny(1),ny(3)) * DELTA(nz(1),nz(3)) )

 if (Qrp.gt.0) then
  V2=V2*(1/Qrp)
 end if

 V3 = DELTA(sz(1),sz(4)) * DELTA(sz(2),sz(1)) * (1 - DELTA(nx(1),nx(4)) * DELTA(ny(1),ny(4)) * DELTA(nz(1),nz(4)) )

 if (Qsp.gt.0) then
  V3=V3*(1/Qsp)
 end if

 V = V1*(V2-V3)

 END FUNCTION MATRIXELEMENT_V

