MODULE interaction_blocks
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



program test 
  use interaction_blocks 
  
  implicit none 
  
  integer :: a,b,i,j,c,d, alpha, Nalpha, Npp, Nhh, dim1, bra,ket, tot_orbs, below_ef, dim2,ket2
  integer, allocatable :: lookup_pp(:,:), lookup_hh(:,:)
  real*8 :: denom, mbpt2, sum1
  
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
  

end program test




