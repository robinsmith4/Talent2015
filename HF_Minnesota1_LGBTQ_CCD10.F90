MODULE CONSTANTS

 Double precision, parameter, public :: pi=3.14159274
 Double precision, parameter, public :: hbar=197.326 !MeV c
 Double precision, parameter, public :: me=0.511 !MeV / c^2
 Double precision, parameter, public :: mn=939.565 !MeV / c^2
 Integer, parameter, public :: NumPart=14

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
	Double precision, external			     :: minnesota_ass
	Double precision, external			     :: FOCK_ME
	Integer, dimension(:), allocatable	   	     :: nx, ny, nz, spin
	Integer						     :: nmax,max_orbits,x,y,z,s,i,j,k,nqs,qnums,t,o,tmin,tmax,a
	Integer, dimension(:,:), allocatable	             :: states
	Double precision, dimension(:), allocatable	     :: E
	Integer, dimension(5)				     :: HOLD
        Double precision				     :: Density
	Double precision				     :: L,Vol,r0
        Double precision				     :: ee,a0,prefac,Kf
	Integer						     :: gs
	Double precision, dimension(:,:), allocatable	     :: H_hf
    	Double precision, dimension(:,:), allocatable	     :: C_hf
	Double precision			   	     :: add
        Integer						     :: m,n,gam,delt
	Double precision				     :: Kin,Pot
	Double precision				     :: Tinf
        Integer						     :: NumPoints,NumPointsC
	Double precision				     :: MBPT_E2
	Integer						     :: NumStates,b,c
	Double precision, dimension(:), allocatable	     :: FockDiag
	Double precision, dimension(:,:,:,:), allocatable    :: T_mat1,T_mat2
	Double precision				     :: Ec1,Ec2
	Integer,dimension(:,:), allocatable		     :: Chan_Indx_dum
	Integer,dimension(:,:), allocatable		     :: Chan_Indx
	Integer,dimension(:), allocatable		     :: Chan_States_Num
        Integer,dimension(:,:,:), allocatable		     :: Chan_States
        Integer						     :: counter,NumChans
        Integer						     :: Nxsum,Nysum,Nzsum,Szsum,Tzsum,Max_ab

 NumPoints = 20
 NumPointsC=20

 tmin=1
 tmax=1

 nmax=1
 max_orbits = 1300
 qnums = 5
! A=14
 gs=2
 Kf=1.0
 r0=1.0

 ee=1.4399764
 a0=52917.72

 prefac=(hbar**2)/(2.*mn)

  allocate ( nx(max_orbits), ny(max_orbits), nz(max_orbits), spin(max_orbits) )
  allocate ( E(max_orbits) )
  allocate ( states(qnums,max_orbits) )
  allocate ( H_hf(max_orbits,max_orbits) )
  allocate ( C_hf(NumPart,max_orbits) )

! do NumPointsC=1,NumPoints

 Density=0.08*(1.0*NumPointsC/NumPoints) !fm^-3

! Density=0.16 !fm^-3

 Vol = (NumPart/Density)

 L = Vol**(1./3.)

 Kf=( (6.*(pi**2))*Density/gs)**(1./3.)

! write(6,*)'Kf= ',Kf

 Tinf=(3*(hbar**2)*(Kf**2))/(10.*mn)

! L=((6*(pi**2)*A)/(gs*Kf**3))**(1./3.)

! write(6,*)L, pi

!----------------------------------------------------------------
!----------------------------------------------------------------

!** SETTING UP THE BASIS **

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
      do t =tmin,tmax


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

 NumStates=i

 allocate ( FockDiag(NumStates) )


 DO j=1,i
  write(6,*)j,states(1,j),states(2,j),states(3,j),states(4,j),states(5,j),E(j)
 END DO

!----------------------------------------------------------------
!----------------------------------------------------------------

!** FOCK optimising (Unused) **



! write(6,*)MATRIXELEMENT_V (1,12,7,10,L,states,max_orbits,qnums)

! First guess for C_hf
! do m=1,A
!  do n=1,i
!   C_hf(m,n)=0
!   if(m.eq.n) C_hf(m,n)=1
!  end do
!  write(6,*)C_hf(m,1),C_hf(m,2),C_hf(m,3),C_hf(m,4),C_hf(m,5)
! end do

! do m=1,i
!  do n=1,i
!   add=0
!   H_hf(m,n) = E(m)*DELTA(m,n)
!   do j=1,A
!    do gam=1,i
!     do delt=1,i
!      add=add+(C_hf(j,gam)*C_hf(j,delt)*MATRIXELEMENT_V(m,gam,n,delt,L,states,max_orbits,qnums))
!      if (MATRIXELEMENT_V(m,gam,n,delt,L,states,max_orbits,qnums).gt.0.00001) then
!       write(6,*)m,gam,n,delt,MATRIXELEMENT_V(m,gam,n,delt,L,states,max_orbits,qnums)
!      endif
!     end do
!    end do
!   end do
!   H_hf(m,n)=H_hf(m,n)+add
!   write(6,*)'Element added ',m,n,H_hf(m,n)
!  end do
! end do

!----------------------------------------------------------------
!----------------------------------------------------------------
!** Calculating The hard tree-fuck energy **!

 Kin=0.0
 Pot=0.0
 do m=1,NumPart
  Kin=Kin+E(m)
  do j=1,NumPart
!   Pot=Pot+0.5*(vmom_minnesota(m,j,m,j,L,states,max_orbits,qnums) &
!      -vmom_minnesota(m,j,j,m,L,states,max_orbits,qnums))
   Pot=Pot+0.5*minnesota_ass(m,j,m,j,L,states,max_orbits,qnums)
  enddo
 enddo

!write(6,*)'Kinetic Energy=',Kin/A,' MeV'
!write(6,*)'Potential energy=',Pot/A,' MeV'
!write(6,*)'Hard tree-fock energy ',(Kin+Pot)/A,' MeV'
!write(6,*)'Inf Kinetic Energy=',Tinf,' MeV'
!write(6,*)'1-T/Tinf ', 1.0 - Tinf/(Kin/A),' MeV'
 
! write(6,*)'Fock Energy'
! write(6,*)Density,(Kin+Pot)/A
! write(6,*)(Kin+Pot)/A



!----------------------------------------------------------------
!----------------------------------------------------------------
!** Calculating the MBPT energy **!

 MBPT_E2=0.0

 do i=1,NumStates
  FockDiag(i)=FOCK_ME(i,i,L,states,max_orbits,qnums,E)
 end do

 do i=1,NumPart
  do j=1,NumPart
   do b=NumPart+1,NumStates
    do c=NumPart+1,NumStates

!     write(6,*)vmom_minnesota(i,j,b,c,L,states,max_orbits,qnums)
     MBPT_E2=MBPT_E2+ ( minnesota_ass(i,j,b,c,L,states,max_orbits,qnums) &
     * minnesota_ass(b,c,i,j,L,states,max_orbits,qnums)  &
     / (FockDiag(i) + FockDiag(j) - FockDiag(b) - FockDiag(c)) )
    end do
   end do
  end do
 end do

 MBPT_E2=0.25*MBPT_E2

! write(6,*)'MBPT Energy'
! write(6,*)Kin,Density,(Kin+MBPT_E2)/A
 write(6,*)Density,(Kin+Pot+MBPT_E2)/NumPart,MBPT_E2

! end do

!----------------------------------------------------------------
!----------------------------------------------------------------
!** Setting up the Channels for CCD **!

 NumChans=((4.*nmax + 1)**3)*3.0*(tmax-tmin+1)

 allocate ( Chan_Indx(NumChans,5) )
 Chan_Indx=0
 counter=0
 do m=-2*nmax,2*nmax
  do i=-2*nmax,2*nmax
   do j=-2*nmax,2*nmax
    do b=-2,2,2
     do c=2*tmin,2*tmax,2
      counter=counter+1
      Chan_Indx(counter,1)=m
      Chan_Indx(counter,2)=i
      Chan_Indx(counter,3)=j
      Chan_Indx(counter,4)=b
      Chan_Indx(counter,5)=c
      write(6,*)counter,m,i,j,b,c
     end do
    end do
   end do
  end do
 end do

 allocate ( Chan_States_Num(NumChans) )
 Chan_States_Num=0
 do a=1,NumStates
  do b=1,NumStates
   Nxsum=states(a,1)+states(b,1)
   Nysum=states(a,2)+states(b,2)
   Nzsum=states(a,3)+states(b,3)
   Szsum=states(a,4)+states(b,4)
   Tzsum=states(a,5)+states(b,5)
   do i=1,NumChans
    if((Nxsum.eq.Chan_Indx(i,1)).and.(Nysum.eq.Chan_Indx(i,2)).and. &
    (Nzsum.eq.Chan_Indx(i,3)).and.(Szsum.eq.Chan_Indx(i,4)).and.(Tzsum.eq.Chan_Indx(i,5))) then
     Chan_States_Num(i)=Chan_States_Num(i)+1
    endif
   end do
  end do
 end do
 
 counter=0
 do i=1,NumChans
  write(6,*)Chan_States_Num(i)
 end do
! write(6,*)counter

 Max_ab=MaxVal(Chan_States_Num)

! write(6,*)Max_ab

 allocate ( Chan_States(NumChans,Max_ab,Max_ab) )


!----------------------------------------------------------------
!----------------------------------------------------------------
!** Calculating the CCD Brueckner HF energy **!

 allocate ( T_mat1(NumPart,NumPart,NumStates,NumStates) )
 allocate ( T_mat2(NumPart,NumPart,NumStates,NumStates) )
 T_mat1=0.0
 T_mat2=0.0

 do i=1,NumPart
  do j=1,NumPart
   do b=NumPart+1,NumStates
    do c=NumPart+1,NumStates
     T_mat1(i,j,b,c)=minnesota_ass(i,j,b,c,L,states,max_orbits,qnums) &
     / (FockDiag(i) + FockDiag(j) - FockDiag(b) - FockDiag(c))
     T_mat2(i,j,b,c)=minnesota_ass(i,j,b,c,L,states,max_orbits,qnums) &
     / (FockDiag(i) + FockDiag(j) - FockDiag(b) - FockDiag(c))
!     write(6,*)T_mat(i,j,b,c)
    end do
   end do
  end do
 end do

 Ec1=0
 Ec2=0

 do i=1,NumPart
  do j=1,NumPart
   do b=NumPart+1,NumStates
    do c=NumPart+1,NumStates
     Ec1=Ec1+( T_mat1(i,j,b,c)*minnesota_ass(i,j,b,c,L,states,max_orbits,qnums) )
    end do
   end do
  end do
 end do

 Ec1=Ec1*0.25

 write(6,*)'Initial E_corr ',Ec1

 do m=1,1000
!  write(6,*)'Starting Iteration ',m
  CALL CCD_Opt(L,FockDiag,NumStates,states,max_orbits,qnums,E,T_mat1,T_mat2)
  do i=1,NumPart
   do j=1,NumPart
    do b=NumPart+1,NumStates
     do c=NumPart+1,NumStates
      Ec2=Ec2+( T_mat2(i,j,b,c)*minnesota_ass(i,j,b,c,L,states,max_orbits,qnums) )
     end do
    end do
   end do
  end do
  T_mat1=T_mat2
  write(6,*)'E_corr iteration ',m,Ec2
  if(abs(Ec2-Ec1).lt.1.0e-4) go to 100
  Ec1=Ec2
  Ec2=0.0
 end do
100 continue
 write(6,*)'Converged E_corr',Ec2



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

    V = vdir/(L**3)
    

  end function vmom_minnesota

!-------------------------------------------------------------------------

 FUNCTION minnesota_ass(p,q,r,s,L,states,max_orbits,qnums) RESULT(Ass)
 USE CONSTANTS

 integer p,q,r,s,max_orbits,qnums
 Double precision L,DELTA,vmom_minnesota
 Double precision Ass
 integer states(qnums,max_orbits)

 Ass=vmom_minnesota(p,q,r,s,L,states,max_orbits,qnums) &
    -vmom_minnesota(p,q,s,r,L,states,max_orbits,qnums)
 
 end function minnesota_ass
!-------------------------------------------------------------------------

 FUNCTION FOCK_ME(p,q,L,states,max_orbits,qnums,E) RESULT(Fpq)
 USE CONSTANTS

 integer p,q,max_orbits,qnums
 Double precision L,DELTA,minnesota_ass
 Double precision Fpq,Tpq,Vpq
 integer states(qnums,max_orbits)
 Double precision E(max_orbits)

 Tpq=E(p)*DELTA(p,q)
 Vpq=0.0

 do i=1,NumPart
  Vpq=Vpq+minnesota_ass(p,i,q,i,L,states,max_orbits,qnums)
 end do

 Fpq=Tpq+Vpq

 end function FOCK_ME
!-------------------------------------------------------------------------
 SUBROUTINE CCD_Opt(L,FockDiag,NumStates,states,max_orbits,qnums,E,T1,T2)
 USE CONSTANTS
 
 integer max_orbits,qnums,i,j,a,b,c,d
 Double precision L,DELTA,minnesota_ass,FOCK_ME
 integer states(qnums,max_orbits)
 Double precision E(max_orbits),FockDiag(NumStates)
 Double precision T1(NumPart,NumPart,NumStates,NumStates),T2(NumPart,NumPart,NumStates,NumStates)

 T_mat2=0.0

! write(6,*)'In CCD opt'

 do i=1,NumPart
  do j=1,NumPart
   do a=NumPart+1,NumStates
    do b=NumPart+1,NumStates
!     write(6,*)'Calculating T_',i,j,a,b
     T2(i,j,a,b)=minnesota_ass(a,b,i,j,L,states,max_orbits,qnums)
     do c=NumPart+1,NumStates
      do d=NumPart+1,NumStates
       T2(i,j,a,b)=T2(i,j,a,b) &
       + 0.5*(minnesota_ass(a,b,c,d,L,states,max_orbits,qnums)*T1(i,j,c,d))
      end do
     end do
     T2(i,j,a,b)=(1./(FockDiag(i)+FockDiag(j)-FockDiag(a)-FockDiag(b)))*T2(i,j,a,b)
    end do
   end do
  end do
 end do

 return
 end
