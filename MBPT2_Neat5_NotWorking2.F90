!****************************************************************************!
!******************************** MODULES ***********************************!

MODULE CONSTANTS

 Double precision, parameter, public 		:: pi=3.14159274
 Double precision, parameter, public 		:: hbar=197.326 !MeV c
 Double precision, parameter, public 		:: me=0.511 !MeV / c^2
 Double precision, parameter, public 		:: mn=939.565 !MeV / c^2
 Double precision, parameter, public 		:: mp=938.272 !MeV / c^2
 Integer, parameter, public 			:: NumPart=14
 Integer, dimension(:,:), allocatable, public 	:: states
 Double precision, allocatable, public 		:: L
 Integer, dimension(:,:), allocatable, public	:: Chan_Indx

END MODULE CONSTANTS

MODULE INTERACTION_BLOCKS
 TYPE, PUBLIC :: block_storage
  REAL*8, DIMENSION(:,:), ALLOCATABLE :: val
 END TYPE block_storage

 TYPE (block_storage), ALLOCATABLE, PUBLIC :: vnn_hhhh(:)
 TYPE (block_storage), ALLOCATABLE, PUBLIC :: vnn_hhhp(:)
 TYPE (block_storage), ALLOCATABLE, PUBLIC :: vnn_hhpp(:)
 TYPE (block_storage), ALLOCATABLE, PUBLIC :: vnn_hphp(:)
 TYPE (block_storage), ALLOCATABLE, PUBLIC :: vnn_hppp(:)
 TYPE (block_storage), ALLOCATABLE, PUBLIC :: vnn_pppp(:)

 TYPE (block_storage), ALLOCATABLE, PUBLIC :: t2_ccm(:)
 TYPE (block_storage), ALLOCATABLE, PUBLIC :: t2_denom(:)
 TYPE (block_storage), ALLOCATABLE, PUBLIC :: t2_ccm_eqn(:)
  
END MODULE INTERACTION_BLOCKS

!****************************************************************************!
!***************************** START PROGRAM ********************************!

PROGRAM BETTERSTORAGE 
 USE INTERACTION_BLOCKS
 USE CONSTANTS
  
 IMPLICIT NONE

!------------------------------ FUNCTIONS ----------------------------------!
 Integer, external 				:: CHANNELCHECK
 Double precision, external 			:: DELTA
 Double precision, external			:: vmom_minnesota
 Double precision, external			:: chp_sigma_dot_sigma
 Double precision, external			:: chp_tau_dot_tau
 Double precision, external			:: minnesota_ass
 Double precision, external			:: FOCK_ME
  
 integer :: a,b,i,j,c,d, alpha, Nalpha, Npp, Nhh, dim1, bra,ket, tot_orbs, below_ef, dim2,ket2
 integer, allocatable :: lookup_pp(:,:), lookup_hh(:,:)
 real*8 :: denom, mbpt2, sum1

!---------------------------- SYSTEM CONSTANTS ------------------------------!

 Integer 		:: tmax,tmin,smax,smin,nmax,max_orbits,qnums,gs !orbits/basis
 Double precision	:: prefac,Density,Kf,Tinf,Vol ! Misc
 Double precision, dimension(:), allocatable ::E
 Integer		:: counter,NumStates,nqs
 Double precision, dimension(:), allocatable :: FockDiag
 Double precision	:: CCD_E

!--------------------------- CHANNEL CONSTANTS ------------------------------!

 Integer,dimension(:), allocatable		:: Chan_States_Num

!****************************************************************************!
!-------------------------- INITIALISE CONSTANTS ----------------------------!

 nmax=1
 tmax=1
 tmin=1
 smin=-1
 smax=1

 qnums=5

 gs=2 !for pnm (change to gs=4 for snm)
 
 if(tmax.ne.tmin) max_orbits=((2*nmax+1)**3)*2*2
 if(tmax.eq.tmin) max_orbits=((2*nmax+1)**3)*2

 prefac=(hbar**2)/(2.*mn)

 Density= 0.12 ! 0.08 !fm^-3

 Vol = (NumPart/Density)

 L = Vol**(1./3.)

 Kf=( (6.*(pi**2))*Density/gs)**(1./3.)

 Tinf=(3*(hbar**2)*(Kf**2))/(10.*mn)

!****************************************************************************!
!---------------------------- INITIALISE BASIS ------------------------------!

 allocate ( E(max_orbits) )
 allocate ( states(qnums,max_orbits) )

 counter=0

 do nqs = 0, nmax*nmax*3

  do a=-nmax,nmax
   do b=-nmax,nmax
    do c=-nmax,nmax
     do d =smin,smax,2
      do i =tmin,tmax,2

	 if (a*a+b*b+c*c /=  nqs)  cycle
	 counter = counter+1

	 states(1,counter) = a
	 states(2,counter) = b
	 states(3,counter) = c
	 states(4,counter) = d
	 states(5,counter) = i

         E(counter)=prefac*((2*pi/L)**2) * (a**2 + b**2 + c**2)

      end do
     end do
    end do
   end do
  end do

 end do

 NumStates=counter

!****************************************************************************!
!-------------------------------- CHANNELS ----------------------------------!

 Nalpha=((4.*nmax + 1)**3)*3*(tmax-tmin+1)

 below_ef = NumPart 
 tot_orbs = max_orbits
  
 ! allocate array for bra and ket lookups.
  allocate( lookup_pp(below_ef+1:tot_orbs, below_ef+1:tot_orbs) )
  allocate( lookup_hh(1:below_ef,1:below_ef) )
  ! allocate interactions for all possible channels
  allocate( vnn_hhpp(Nalpha) ) 
  allocate( vnn_pppp(Nalpha) ) 
  allocate( vnn_hhhh(Nalpha) ) 
  allocate( t2_ccm(Nalpha) ) 
  allocate( t2_denom(Nalpha) )
  allocate( t2_ccm_eqn(Nalpha) ) 
    
  lookup_pp = 0
  lookup_hh = 0

!****************************************************************************!
!------------------------ CALCULATE CHANNEL DETAILS -------------------------!

 allocate ( Chan_Indx(Nalpha,qnums) )
 Chan_Indx=0
 counter=0
 do a=-2*nmax,2*nmax
  do b=-2*nmax,2*nmax
   do c=-2*nmax,2*nmax
    do d=2*smin,2*smax,2
     do i=2*tmin,2*tmax,2
      counter=counter+1
      Chan_Indx(counter,1)=a
      Chan_Indx(counter,2)=b
      Chan_Indx(counter,3)=c
      Chan_Indx(counter,4)=d
      Chan_Indx(counter,5)=i
!      write(6,*)counter,a,b,c,d,i
     end do
    end do
   end do
  end do
 end do

 Nalpha=counter

 allocate ( FockDiag(NumStates) )
 do i=1,NumStates
  FockDiag(i)=FOCK_ME(i,i,NumStates,E)
 end do

!****************************************************************************!
!------------------------- CALCULATE CHANNEL SIZES --------------------------!

  ! Before here, we have calculated each allowed 'channel' (sum q nums) and the number of channels
  ! Now to allocate and setup each matrix for each channel

 do alpha = 1, Nalpha
 
 write(6,*)alpha

  ! compute bra dimension of vhhpp (depends on channel alpha)
  dim1 = 0
  do i = 1, below_ef
   do j = 1, below_ef
    ! check if this combo of i and j belong to alpha
    if(ChannelCheck(alpha,i,j).eq.1) then
     dim1 = dim1 + 1
    endif
   end do
  end do


! compute ket dimension of vpppp (depends on channel alpha)
  dim2 = 0
  do i = below_ef+1, tot_orbs 
   do j = below_ef+1, tot_orbs
    ! check if this combo of i and j belong to alpha
    if(ChannelCheck(alpha,i,j).eq.1) then
     dim2 = dim2 + 1
    endif
   end do
  end do

!****************************************************************************!
!----------------- ALLOCATE SPACE FOR VNN AND T2 ETC... ---------------------!

  allocate( vnn_pppp(alpha)%val(dim2, dim2) )
  allocate( vnn_hhpp(alpha)%val(dim1, dim2) )
  allocate( vnn_hhhh(alpha)%val(dim1, dim1) )
  allocate( t2_ccm(alpha)%val(dim2, dim1) )
  allocate( t2_denom(alpha)%val(dim2, dim1) )
  allocate( t2_ccm_eqn(alpha)%val(dim2, dim1) )

  vnn_pppp(alpha)%val = 0.d0
  vnn_hhpp(alpha)%val = 0.d0
  vnn_hhhh(alpha)%val = 0.d0
  t2_ccm(alpha)%val = 0.d0
  t2_ccm_eqn(alpha)%val = 0.d0

!****************************************************************************!
!------------------------ CALCULATE MATRIX ELEMENTS -------------------------!
! compute pppp matrix elements

  bra=0
  ket=0
  do a = below_ef+1,tot_orbs 
   do b = below_ef+1,tot_orbs
    if(ChannelCheck(alpha,a,b).eq.1) then
     bra=bra+1
     ket=0
     do i = below_ef+1,tot_orbs 
      do j = below_ef+1,tot_orbs
       if(ChannelCheck(alpha,i,j).eq.1) then
        ket = ket + 1
        vnn_pppp(alpha)%val(bra,ket) = minnesota_ass(a,b,i,j) ! Matrix element a,b,i,j
       end if
      end do
     end do
    end if
   end do
  end do

! compute hhhh matrix elements

  bra=0
  ket=0
  do a = 1,below_ef
   do b = 1,below_ef
    if(ChannelCheck(alpha,a,b).eq.1) then
     bra=bra+1
     ket=0
     do i = 1,below_ef
      do j = 1,below_ef
       if(ChannelCheck(alpha,i,j).eq.1) then
        ket = ket + 1
        vnn_hhhh(alpha)%val(bra,ket) = minnesota_ass(a,b,i,j) ! Matrix element a,b,i,j
       end if
      end do
     end do
    end if
   end do
  end do

! compute hhpp matrix elements       
  bra = 0
  ket = 0 
  do a = below_ef+1,tot_orbs 
   do b = below_ef+1,tot_orbs
    if(ChannelCheck(alpha,a,b).eq.1) then
     bra = bra + 1 
     lookup_pp(a,b) = bra
     ket = 0 
     do i = 1, below_ef
      do j = 1, below_ef
       if(ChannelCheck(alpha,i,j).eq.1) then
        ket = ket + 1 
        lookup_hh(i,j) = ket
        vnn_hhpp(alpha)%val(ket,bra) = minnesota_ass(i,j,a,b) !bra*ket 
        t2_ccm(alpha)%val(bra,ket) = minnesota_ass(a,b,i,j) &
        /(FockDiag(i) + FockDiag(j) - FockDiag(a) - FockDiag(b)) !bra*ket*2.
        t2_denom(alpha)%val(bra,ket) = &
        1.0/(FockDiag(i) + FockDiag(j) - FockDiag(a) - FockDiag(b))
        !t2_denom(alpha)%val(bra,ket) = 1
        t2_ccm_eqn(alpha)%val(bra,ket) = 0.0
       end if
      end do
     end do
    end if
   end do
  end do


 end do

!****************************************************************************!
!----------------------------------- MBPT -----------------------------------!
 ! compute MBPT2
 sum1 = 0.d0 
 do alpha = 1, Nalpha
  do bra = 1, size( vnn_hhpp(alpha)%val, 1) ! dim1 and dim2 ?
   do ket = 1, size( vnn_hhpp(alpha)%val, 2)
    sum1=sum1+(vnn_hhpp(alpha)%val(bra,ket)*t2_ccm(alpha)%val(ket,bra))
   end do
  end do
 end do
 sum1=sum1/4.0

 write(6,*)'MBPT2 =',sum1

!****************************************************************************!
!----------------------------------- CCD -----------------------------------!

 do alpha = 1, Nalpha

 ! compute diagram Vpphh
  sum1 = 0.d0
  do bra = 1, size(t2_ccm(alpha)%val, 1) ! dim1 and dim2 ?
   do ket = 1, size(t2_ccm(alpha)%val, 2)
    t2_ccm_eqn(alpha)%val(bra,ket) = t2_ccm_eqn(alpha)%val(bra,ket) + &
    (vnn_hhpp(alpha)%val(bra,ket)*t2_denom(alpha)%val(bra,ket))
   end do
  end do

 ! compute diagram Vpppp*T2
  do bra = 1, size(t2_ccm(alpha)%val, 1) ! dim1 and dim2 ?
   do ket = 1, size(t2_ccm(alpha)%val, 2)
    sum1 = 0.d0 
    do ket2 = 1, size( vnn_pppp(alpha)%val, 1)
     sum1 = sum1 + 0.5*vnn_pppp(alpha)%val(bra,ket2)*t2_ccm(alpha)%val(ket2,ket) &
     *t2_denom(alpha)%val(bra,ket)
    end do
    t2_ccm_eqn(alpha)%val(bra,ket) = t2_ccm_eqn(alpha)%val(bra,ket)+sum1
   end do
  end do

 ! compute diagram Vhhhh*T2

  do bra = 1, size(t2_ccm(alpha)%val, 1) ! dim1 and dim2 ?
   do ket = 1, size(t2_ccm(alpha)%val, 2)
    sum1 = 0.d0 
    do ket2 = 1, size( vnn_hhhh(alpha)%val, 1)
     sum1 = sum1 + 0.5*vnn_hhhh(alpha)%val(bra,ket2)*t2_ccm(alpha)%val(ket2,ket) &
     *t2_denom(alpha)%val(bra,ket)
    end do
    t2_ccm_eqn(alpha)%val(bra,ket)=t2_ccm_eqn(alpha)%val(bra,ket) + sum1
   end do
  end do

 end do

  ! example ccd energy
 CCD_E = 0.0
 do alpha = 1, Nalpha 
  do bra = 1, size(vnn_hhpp(alpha)%val, 1)
   do ket = 1, size(vnn_hhpp(alpha)%val, 2)
    CCD_E = CCD_E+0.25*vnn_hhpp(alpha)%val(bra,ket)*t2_ccm(alpha)%val(ket,bra)     
   end do
  end do
 end do
 write(6,*) 'Initial CCD correlation energy', CCD_E

   ! example ccd energy
 CCD_E = 0.0
 do alpha = 1, Nalpha 
  do bra = 1, size(vnn_hhpp(alpha)%val, 1)
   do ket = 1, size(vnn_hhpp(alpha)%val, 2)
    CCD_E = CCD_E+0.25*vnn_hhpp(alpha)%val(bra,ket)*t2_ccm_eqn(alpha)%val(ket,bra)     
   end do
  end do
 end do
 write(6,*) 'CCD correlation energy', CCD_E

  
100 continue
END PROGRAM BETTERSTORAGE

!-------------------------------------------------------------------------

 FUNCTION CHANNELCHECK(alpha,i,j) RESULT(check)
  USE CONSTANTS
!  USE INTERACTION_BLOCKS

  integer alpha,i,j,check
  integer Nxsum,Nysum,Nzsum,Szsum,Tzsum

  Nxsum=states(1,i)+states(1,j)
  Nysum=states(2,i)+states(2,j)
  Nzsum=states(3,i)+states(3,j)
  Szsum=states(4,i)+states(4,j)
  Tzsum=states(5,i)+states(5,j)

  check=0

  if((Nxsum.eq.Chan_Indx(alpha,1)).and.(Nysum.eq.Chan_Indx(alpha,2)).and. &
     (Nzsum.eq.Chan_Indx(alpha,3)).and.(Szsum.eq.Chan_Indx(alpha,4)).and. &
     (Tzsum.eq.Chan_Indx(alpha,5))) then
   check=1
  endif

 end function CHANNELCHECK
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

 FUNCTION vmom_minnesota(p,q,r,s) RESULT(V)
 USE CONSTANTS
  
    implicit none
!    integer :: max_orbits,qnums
!    integer :: states(qnums,max_orbits)
    Double precision :: V
    Double precision :: DELTA,chp_sigma_dot_sigma,chp_tau_dot_tau
    INTEGER :: p,q,r,s, m1,m2,m3,m4, spin, iph, t1,t2,t3,t4, Tiso
    INTEGER :: nx1, ny1, nz1, nx2, ny2, nz2, nx3, ny3, nz3, nx4, ny4, nz4 
    REAL*8 :: k1(3), k2(3), k3(3), k4(3), kmean(3) 
    REAL*8 :: qtrans(3), prel(3), pprel(3)
    REAL*8 :: q2, p2, qabs, pp2, vdir, vexc, cg1, cg2 
    REAL*8 :: nucleon_mass, relativity_factor !delta,
    REAL*8 :: v0r, v0t, v0s, kr,kt,ks, vr,vt,vs,vs_ex,vcentral, vsigma, spin_exc1, spin_exc2
    
    V=0.0

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

 FUNCTION minnesota_ass(p,q,r,s) RESULT(Ass)
 USE CONSTANTS

 integer p,q,r,s
 Double precision DELTA,vmom_minnesota
 Double precision Ass
! integer states(qnums,max_orbits)

 Ass=vmom_minnesota(p,q,r,s) &
    -vmom_minnesota(p,q,s,r)
 
 end function minnesota_ass
!-------------------------------------------------------------------------

 FUNCTION FOCK_ME(p,q,NumStates,E) RESULT(Fpq)
 USE CONSTANTS

 integer p,q,NumStates
 Double precision DELTA,minnesota_ass
 Double precision Fpq,Tpq,Vpq
 Double precision E(NumStates)

 Tpq=E(p)*DELTA(p,q)
 Vpq=0.0

 do i=1,NumPart
  Vpq=Vpq+minnesota_ass(p,i,q,i)
 end do

 Fpq=Tpq+Vpq

 end function FOCK_ME
!-------------------------------------------------------------------------

