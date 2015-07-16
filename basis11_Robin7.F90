MODULE CONSTANTS

 Double precision, parameter, public :: pi=3.14159274
 Double precision, parameter, public :: hbar=197.326 !MeV c
 Double precision, parameter, public :: me=0.511 !MeV / c^2

END MODULE CONSTANTS

PROGRAM BASIS
USE CONSTANTS

 IMPLICIT NONE

        Double precision, external 			     :: DELTA
	Double precision, external 			     :: DELTA2
        Double precision, external			     :: MATRIXELEMENT_V
	Double precision, external			     :: vmom_minnesota
	Double precision, external			     :: chp_sigma_dot_sigma
	Double precision, external			     :: chp_tau_dot_tau
        Double precision, external			     :: QTRANS
	Integer, dimension(:), allocatable	   	     :: nx, ny, nz, spin
	Integer						     :: nmax,max_orbits,x,y,z,s,i,j,k,nqs,qnums,t,o
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
      do t =-1,1,2


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
	 states(5,i) = t

         E(i)=prefac*((2*pi/L)**2) * (nx(i)**2 + ny(i)**2 + nz(i)**2)

      end do
     end do
    end do
   end do
  end do

 end do


! DO j=1,i
!  write(6,*)j,states(1,j),states(2,j),states(3,j),states(4,j),states(5,j),E(j)
! END DO

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

! do m=-1,1,2
!  do n=-1,1,2
!   do gam=-1,1,2
!    do o=-1,1,2
!     write(6,*)m,n,gam,o,chp_tau_dot_tau_mtx(m,n,gam,o)
!    end do
!   end do
!  end do
! end do

 do j=1,14
  do m=1,14
    write(6,*)vmom_minnesota(m,j,m,j,L,states,max_orbits,qnums)
!   write(6,*)MATRIXELEMENT_V(m,j,m,j,L,states,max_orbits,qnums)
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

!-------------------------------------------------------------------------

 FUNCTION chp_sigma_dot_sigma(ms1,ms2,ms3,ms4) RESULT(res)
    
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    REAL*8 :: res, delta
    complex*16 :: res1
    COMPLEX*16 :: chi1(2), chi2(2), chi3(2), chi4(2)
    INTEGER :: i1 
    COMPLEX*16 :: sigma_x(2,2),sigma_y(2,2),sigma_z(2,2)

    sigma_x(1,1) = (0.0, 0.0)
    sigma_x(1,2) = (1.0, 0.0)
    sigma_x(2,1) = (1.0, 0.0)
    sigma_x(2,2) = (0.0, 0.0)

    sigma_y(1,1) = (0.0, 0.0)
    sigma_y(1,2) = (0.0, -1.0)
    sigma_y(2,1) = (0.0, 1.0)
    sigma_y(2,2) = (0.0, 0.0)

    sigma_z(1,1) = (1.0, 0.0)
    sigma_z(1,2) = (0.0, 0.0)
    sigma_z(2,1) = (0.0, 0.0)
    sigma_z(2,2) = (-1.0, 0.0)
    
    res = 0.0D0
    chi1 = 0.d0; chi2 = 0.d0; chi3 = 0.d0; chi4 = 0.d0
    i1 = nint(1.5-0.5*ms1) 
    chi1(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms2) 
    chi2(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms3) 
    chi3(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms4) 
    chi4(i1) = 1.d0 
    
    
    res1 = dot_product(chi1, matmul( sigma_x, chi3)) * dot_product(chi2, matmul( sigma_x, chi4)) &
         + dot_product(chi1, matmul( sigma_y, chi3)) * dot_product(chi2, matmul( sigma_y, chi4)) &
         + dot_product(chi1, matmul( sigma_z, chi3)) * dot_product(chi2, matmul( sigma_z, chi4)) 
    res = REAL(res1)
    
    
 END FUNCTION chp_sigma_dot_sigma

!-------------------------------------------------------------------------

 FUNCTION chp_tau_dot_tau(ms1,ms2,ms3,ms4) RESULT(res)
    
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    REAL*8 :: res, delta
    complex*16 :: res1
    COMPLEX*16 :: chi1(2), chi2(2), chi3(2), chi4(2)
    INTEGER :: i1 
    COMPLEX*16 :: sigma_x(2,2),sigma_y(2,2),sigma_z(2,2)

    sigma_x(1,1) = (0.0, 0.0)
    sigma_x(1,2) = (1.0, 0.0)
    sigma_x(2,1) = (1.0, 0.0)
    sigma_x(2,2) = (0.0, 0.0)

    sigma_y(1,1) = (0.0, 0.0)
    sigma_y(1,2) = (0.0, -1.0)
    sigma_y(2,1) = (0.0, 1.0)
    sigma_y(2,2) = (0.0, 0.0)

    sigma_z(1,1) = (1.0, 0.0)
    sigma_z(1,2) = (0.0, 0.0)
    sigma_z(2,1) = (0.0, 0.0)
    sigma_z(2,2) = (-1.0, 0.0)
    
    res = 0.0D0
    chi1 = 0.d0; chi2 = 0.d0; chi3 = 0.d0; chi4 = 0.d0
    i1 = nint(1.5-0.5*ms1) 
    chi1(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms2) 
    chi2(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms3) 
    chi3(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms4) 
    chi4(i1) = 1.d0 
    
    
    res1 = dot_product(chi1, matmul( sigma_x, chi3)) * dot_product(chi2, matmul( sigma_x, chi4)) &
         + dot_product(chi1, matmul( sigma_y, chi3)) * dot_product(chi2, matmul( sigma_y, chi4)) &
         + dot_product(chi1, matmul( sigma_z, chi3)) * dot_product(chi2, matmul( sigma_z, chi4)) 
    res = REAL(res1)
    
    
 END FUNCTION chp_tau_dot_tau

!-------------------------------------------------------------------------

 FUNCTION vmom_minnesota(p,q,r,s,L,states,max_orbits,qnums) RESULT(V)
  USE CONSTANTS
!    USE single_particle_orbits
!    USE constants
!    use chiral_constants
  
    implicit none
    integer :: max_orbits,qnums
    integer :: states(qnums,max_orbits)
    Double precision :: L,V
    Double precision :: DELTA,chp_sigma_dot_sigma,chp_tau_dot_tau
    INTEGER :: p,q,r,s, m1,m2,m3,m4, spin, iph, t1,t2,t3,t4, Tiso
    INTEGER :: nx1, ny1, nz1, nx2, ny2, nz2, nx3, ny3, nz3, nx4, ny4, nz4 
    REAL*8 :: k1(3), k2(3), k3(3), k4(3), kmean(3) 
    REAL*8 :: qtrans(3), prel(3), pprel(3)
    REAL*8 :: q2, p2, qabs, pp2, vdir, vexc, cg1, cg2 
    REAL*8 :: nucleon_mass, relativity_factor !delta,
    REAL*8 :: v0r, v0t, v0s, kr,kt,ks, vr,vt,vs,vs_ex,vcentral, vsigma, spin_exc1, spin_exc2
  
    nx1 = states(1,p)
    ny1 = states(2,p)
    nz1 = states(3,p)
    nx2 = states(1,q)
    ny2 = states(2,q)
    nz2 = states(3,q)
    nx3 = states(1,r)
    ny3 = states(2,r)
    nz3 = states(3,r)
    nx4 = states(1,s)
    ny4 = states(2,s)
    nz4 = states(3,s)
    ! 
    ! Conservation of linear momentum
    !
    if ( nx1 + nx2 /= nx3 + nx4 ) return 
    if ( ny1 + ny2 /= ny3 + ny4 ) return 
    if ( nz1 + nz2 /= nz3 + nz4 ) return 
  
    k1(1) = (2.*pi/L)*nx1
    k1(2) = (2.*pi/L)*ny1
    k1(3) = (2.*pi/L)*nz1
    k2(1) = (2.*pi/L)*nx2
    k2(2) = (2.*pi/L)*ny2
    k2(3) = (2.*pi/L)*nz2
    k3(1) = (2.*pi/L)*nx3
    k3(2) = (2.*pi/L)*ny3
    k3(3) = (2.*pi/L)*nz3
    k4(1) = (2.*pi/L)*nx4
    k4(2) = (2.*pi/L)*ny4
    k4(3) = (2.*pi/L)*nz4
  
    ! 
    ! conservation of spin and isospin 
    !
  
    m1 = states(4,p)
    m2 = states(4,q)
    m3 = states(4,r)
    m4 = states(4,s)

    t1 = states(5,p)
    t2 = states(5,q)
    t3 = states(5,r) 
    t4 = states(5,s)

    if ( m1 + m2 /= m3 + m4 ) return
    if ( t1 + t2 /= t3 + t4 ) return
    
    ! 
    ! RELATIVE MOMENTA <prel |v| pprel > 
    ! 
    prel  = 0.5d0*(k1-k2)
    pprel = 0.5d0*(k3-k4)

    !
    ! momentum transfer 
    !
    qtrans = prel - pprel
    q2 = sum(qtrans*qtrans) 
    
    v0r=200.0  ! MeV
    v0t=178.0  ! MeV
    v0s=91.85  ! MeV
    kr=1.487  ! fm**-2
    kt=0.639  ! fm**-2
    ks=0.465  ! fm**-2

  
    ! r-space 
    !vr=v0r*exp(-kr*rr**2)
    !vt=-v0t*exp(-kt*rr**2)
    !vs=-v0s*exp(-ks*rr**2)
    
    vr =  v0r * pi**1.5d0 * exp(-q2/(4.d0*kr) )/ (kr**1.5d0) 
    vt = -v0t * pi**1.5d0 * exp(-q2/(4.d0*kt) )/ (kt**1.5d0)
    vs = -v0s * pi**1.5d0 * exp(-q2/(4.d0*ks) )/ (ks**1.5d0)
   
    
    vr = vr * (3.d0*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4) - & 
         chp_tau_dot_tau(t1,t2,t3,t4)*delta(m1,m3)*delta(m2,m4) - & 
         chp_sigma_dot_sigma(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4) - & 
         chp_tau_dot_tau(t1,t2,t3,t4)*chp_sigma_dot_sigma(m1,m2,m3,m4) )/8.d0
    
    vs = vs * (3.d0*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4) + & 
         chp_tau_dot_tau(t1,t2,t3,t4)*delta(m1,m3)*delta(m2,m4) - & 
         3.d0*chp_sigma_dot_sigma(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4) - & 
         chp_tau_dot_tau(t1,t2,t3,t4)*chp_sigma_dot_sigma(m1,m2,m3,m4) )/16.d0
    
    vt = vt * (3.d0*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4) - & 
         3.d0*chp_tau_dot_tau(t1,t2,t3,t4)*delta(m1,m3)*delta(m2,m4) + & 
         chp_sigma_dot_sigma(m1,m2,m3,m4)*delta(t1,t3)*delta(t2,t4) - & 
         chp_tau_dot_tau(t1,t2,t3,t4)*chp_sigma_dot_sigma(m1,m2,m3,m4) )/16.d0
    
   
    vdir = vs+vr+vt
    
!    vmom_minnesota = vdir/volume


    V = vdir/(L**3)
    

  end function vmom_minnesota

!-------------------------------------------------------------------------
