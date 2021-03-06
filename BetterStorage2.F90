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
 TYPE (block_storage), ALLOCATABLE, PUBLIC :: t2_ccm_eqn(:)
  
END MODULE INTERACTION_BLOCKS

!****************************************************************************!
!***************************** START PROGRAM ********************************!

PROGRAM BETTERSTORAGE 
 USE INTERACTION_BLOCKS
 USE CONSTANTS
  
 IMPLICIT NONE
  
 integer :: a,b,i,j,c,d, alpha, Nalpha, Npp, Nhh, dim1, bra,ket, tot_orbs, below_ef, dim2,ket2
 integer, allocatable :: lookup_pp(:,:), lookup_hh(:,:)
 real*8 :: denom, mbpt2, sum1

!****************************************************************************!
!---------------------------- SYSTEM CONSTANTS ------------------------------!

 Integer 		:: tmax,tmin,smax,smin,nmax,max_orbits,qnums,gs !orbits/basis
 Double precision	:: prefac,Density,Kf,Tinf,Vol ! Misc
 Double precision, dimension(:), allocatable ::E
 Integer		:: counter,NumStates,nqs


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

 Density=0.08 !fm^-3

 Vol = (NumPart/Density)

 L = Vol**(1./3.)

 write(6,*)L

 Kf=( (6.*(pi**2))*Density/gs)**(1./3.)

 Tinf=(3*(hbar**2)*(Kf**2))/(10.*mn)

 write(6,*)max_orbits

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

         write(6,*)counter,a,b,c,d,i,E(counter)

      end do
     end do
    end do
   end do
  end do

 end do

 NumStates=i
  
 Nalpha = 2 
  
 below_ef = 2 
 tot_orbs = 8
  
 ! allocate array for bra and ket lookups.
  allocate( lookup_pp(below_ef+1:tot_orbs, below_ef+1:tot_orbs) )
  allocate( lookup_hh(1:below_ef,1:below_ef) )
  ! allocate interactions for all possible channels
  allocate( vnn_hhpp(Nalpha) ) 
  allocate( vnn_pppp(Nalpha) ) 
  allocate( t2_ccm(Nalpha) ) 
  allocate( t2_ccm_eqn(Nalpha) ) 
    
  lookup_pp = 0
  lookup_hh = 0
  ! Before here, need to calculate each allowed 'channel' (sum q nums and the number of channels)
  ! allocate and setup each matrix for each channel
  do alpha = 1, Nalpha 
     
     ! compute bra dimension of vhhpp (depends on channel alpha)
     dim1 = 0
     do i = 1, below_ef
        do j = 1, below_ef
           ! check if this combo of i and j belong to alpha
           dim1 = dim1 + 1 
        end do
     end do

     ! compute ket dimension of vpppp (depends on channel alpha)
     dim2 = 0
     do i = below_ef+1, tot_orbs 
        do j = below_ef+1, tot_orbs
           ! check if this combo of i and j belong to alpha
           dim2 = dim2 + 1 
        end do
     end do
     
     allocate( vnn_pppp(alpha)%val(dim2, dim2) )
     allocate( vnn_hhpp(alpha)%val(dim1, dim2) )
     allocate( t2_ccm(alpha)%val(dim2, dim1) )
     allocate( t2_ccm_eqn(alpha)%val(dim2, dim1) )

     vnn_pppp(alpha)%val = 0.d0
     vnn_hhpp(alpha)%val = 0.d0
     t2_ccm(alpha)%val = 0.d0
     t2_ccm_eqn(alpha)%val = 0.d0

! compute pppp matrix elements     
     do bra = 1, dim2 
        do ket = 1, dim2 
           vnn_pppp(alpha)%val(bra,ket) = bra**2/ket
           ! Need the i,j,a,b here
        end do
     end do
! compute hhpp matrix elements       
     bra = 0 
     ! do bra=1,dim2
     do a = below_ef+1,tot_orbs 
        do b = below_ef+1,tot_orbs
           bra = bra + 1 
           lookup_pp(a,b) = bra
           
           ket = 0 
           do i = 1, below_ef
              do j = 1, below_ef
                 ket = ket + 1 
                 lookup_hh(i,j) = ket

                 write(6,*)i,j,lookup_hh(i,j)
                 
                 vnn_hhpp(alpha)%val(ket,bra) = bra*ket 
                 t2_ccm(alpha)%val(bra,ket) = bra*ket*2. 
                 ! Calculate v_minn using i,j,a,b here
                 
                 
              end do
           end do
        end do
     end do
  end do

  

  ! compute diagram Vpppp*T2 
  do alpha = 1, Nalpha 
     
     
     do bra = 1, size(t2_ccm(alpha)%val, 1) ! dim1 and dim2 ?
        do ket = 1, size(t2_ccm(alpha)%val, 2)
           
           sum1 = 0.d0 
           do ket2 = 1, size( vnn_pppp(alpha)%val, 1)

              sum1 = sum1 + 0.5*vnn_pppp(alpha)%val(bra,ket2)*t2_ccm(alpha)%val(ket2,ket)
           end do
           
           
           t2_ccm_eqn(alpha)%val(bra,ket) = sum1

!!$           t2_ccm_eqn(alpha)%val = matmul(0.5*vnn_pppp(alpha)%val,t2_ccm(alpha)%val)
           
        end do
     end do
  end do

  ! expample ccd energy
  mbpt2 = 0. 
  do alpha = 1, Nalpha 
     
     do bra = 1, size(vnn_hhpp(alpha)%val, 1)
        do ket = 1, size(vnn_hhpp(alpha)%val, 2)
           
           denom = bra
           mbpt2 = mbpt2 + 0.25*vnn_hhpp(alpha)%val(bra,ket)*t2_ccm_eqn(alpha)%val(ket,bra)/denom
                 
        end do
     end do
  end do
  write(6,*) 'CCD correlation energy', mbpt2
  

END PROGRAM BETTERSTORAGE




