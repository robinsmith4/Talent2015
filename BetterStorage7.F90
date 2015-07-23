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
 TYPE (block_storage), ALLOCATABLE, PUBLIC :: t2_ccm_eqn(:)
  
END MODULE INTERACTION_BLOCKS

!****************************************************************************!
!***************************** START PROGRAM ********************************!

PROGRAM BETTERSTORAGE 
 USE INTERACTION_BLOCKS
 USE CONSTANTS
  
 IMPLICIT NONE

!Functions
 Integer, external 			     :: CHANNELCHECK
  
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
!--------------------------- CHANNEL CONSTANTS ------------------------------!

! Integer,dimension(:,:), allocatable		:: Chan_Indx
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

!         write(6,*)counter,a,b,c,d,i,E(counter)

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
  
! Nalpha = 2
 
! below_ef = 2
! tot_orbs = 8
 below_ef = NumPart 
 tot_orbs = max_orbits
  
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
      write(6,*)counter,a,b,c,d,i
     end do
    end do
   end do
  end do
 end do


!****************************************************************************!!
! Nalpha=10

  ! Before here, we have calculated each allowed 'channel' (sum q nums and the number of channels)
  ! Now to allocate and setup each matrix for each channel

  do alpha = 1, Nalpha
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
!     if(dim1.gt.0) write(6,*)'Alpha ',alpha,'Dimension 1 ',dim1

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
!     if(dim2.gt.0) write(6,*)'Alpha ',alpha,'Dimension 2 ',dim2
     
     allocate( vnn_pppp(alpha)%val(dim2, dim2) )
     allocate( vnn_hhpp(alpha)%val(dim1, dim2) )
     allocate( t2_ccm(alpha)%val(dim2, dim1) )
     allocate( t2_ccm_eqn(alpha)%val(dim2, dim1) )

     vnn_pppp(alpha)%val = 0.d0
     vnn_hhpp(alpha)%val = 0.d0
     t2_ccm(alpha)%val = 0.d0
     t2_ccm_eqn(alpha)%val = 0.d0

! **** The work is completed up to here ****
!****************************************************************************!
!------------------------ CALCULATE MATRIX ELEMENTS -------------------------!
! compute pppp matrix elements
     bra=0
     do a = below_ef+1,tot_orbs 
        do b = below_ef+1,tot_orbs
	     if(ChannelCheck(alpha,a,b).eq.1) then
              bra=bra+1
	     end if
             ket=0
             do i = below_ef+1,tot_orbs 
                  do j = below_ef+1,tot_orbs
		     	if(ChannelCheck(alpha,i,j).eq.1) then
              		 ket = ket + 1
			 vnn_pppp(alpha)%val(bra,ket) = bra**2/ket
	     		end if
		end do
           end do
        end do
     end do

!     do bra = 1, dim2 
!        do ket = 1, dim2 
!           vnn_pppp(alpha)%val(bra,ket) = bra**2/ket
           ! Need the i,j,a,b here
!        end do
!B     end do
! compute hhpp matrix elements       
     bra = 0 
     ! do bra=1,dim2
     do a = below_ef+1,tot_orbs 
        do b = below_ef+1,tot_orbs
           if(ChannelCheck(alpha,a,b).eq.1) then
            bra = bra + 1 
            lookup_pp(a,b) = bra
           end if
           
           ket = 0 
           do i = 1, below_ef
              do j = 1, below_ef
!		 write(6,*)a,b,i,j
!                 write(6,*)i,j,lookup_hh(i,j)
                 if(ChannelCheck(alpha,i,j).eq.1) then
		  ket = ket + 1 
                  lookup_hh(i,j) = ket
                  vnn_hhpp(alpha)%val(ket,bra) = bra*ket 
                  t2_ccm(alpha)%val(bra,ket) = bra*ket*2. 
                 ! Calculate v_minn using i,j,a,b here
                 end if
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

  ! example ccd energy
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

  
100 continue
END PROGRAM BETTERSTORAGE

!-------------------------------------------------------------------------

 FUNCTION CHANNELCHECK(alpha,i,j) RESULT(check)
  USE CONSTANTS
!  USE INTERACTION_BLOCKS

  integer alpha,i,j,check
  integer Nxsum,Nysum,Nzsum,Szsum,Tzsum
!  write(6,*)alpha,i,j
  Nxsum=states(1,i)+states(1,j)
  Nysum=states(2,i)+states(2,j)
  Nzsum=states(3,i)+states(3,j)
  Szsum=states(4,i)+states(4,j)
  Tzsum=states(5,i)+states(5,j)
!  write(6,*)Nxsum,Nysum,Nzsum,Szsum,Tzsum,Chan_Indx(alpha,1)

  check=0

  if((Nxsum.eq.Chan_Indx(alpha,1)).and.(Nysum.eq.Chan_Indx(alpha,2)).and. &
     (Nzsum.eq.Chan_Indx(alpha,3)).and.(Szsum.eq.Chan_Indx(alpha,4)).and. &
     (Tzsum.eq.Chan_Indx(alpha,5))) then
   check=1
  endif

 end function CHANNELCHECK
!-------------------------------------------------------------------------


