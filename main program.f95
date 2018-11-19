!****************************************************************************************************
        
!****************************************************************************************************
                                                                                               !main program starts 
!****************************************************************************************************
program cms
use common_variables
implicit none
include 'mpif.h'
integer ista,iend,ierr,iproc,nproc
integer config, complete

        integer :: bin1 = 0
        integer :: mbin1=0
        integer :: bin2 = 0
        integer :: mbin2=0
        integer :: bin3 = 0
        integer :: mbin3=0
        integer :: bin4 = 0
        integer :: mbin4=0
        integer :: bin5 = 0
        integer :: mbin5=0
        integer :: bin6 = 0
        integer :: mbin6=0
        integer :: bin7 = 0
        integer :: mbin7=0
        integer :: bin8 = 0
        integer :: mbin8=0
        integer :: bin9 = 0
        integer :: mbin9=0
        
    real*8  start, finish
    real*8 :: ab_x(spec_step) = 0.D0
    real*8 :: mab_x(spec_step) = 0.D0
    real*8 :: ab_y(spec_step) = 0.D0
    real*8 :: mab_y(spec_step) = 0.D0
    real*8 :: sumoscx = 0.d0, sumoscy = 0.d0
    real*8 :: msumoscx = 0.d0, msumoscy = 0.d0
    real*8, allocatable :: cohfxn(:,:,:), cohtmp(:,:,:)
    real*8, allocatable :: mcohfxn(:,:,:)    
    character*128 sstring
    
  
    !set the parameters (local)
    call set_para()
    
  
   
    
    !do the indexing (exciton_common_local)
    
!****************************************************************************************************
  !We need to first index the lattice and the basis state, that is the one particle and two particle states
    call index_lattice()
    if ( one_state )   call index_1p()
    if ( two_state )   call index_2p()    
!****************************************************************************************************      


 !**************************output the total no. of one particle and two particle basis sets for our check************!
    print*, '***************************'
    print*, 'kount   :', kount
    print*, 'kount 1p:', kount_1p
    print*, 'kount 2p:', kount_2p
    print*, '***************************'
!****************************************************************************************************      
  
    
!****************************************************************************************************
                                                  !We need to set the Frank Condon factors for the vibrational overlap
    call set_fctable()
!****************************************************************************************************        
    
    !output the parameters (local)
    call para_out_cms()
    
 !****************************************************************************************************
                     !Here we are calling the subroutines for the diagonal and off-diagonal disorder(paracrystalline_SS)
    
    if (diagonal) call disorder_elements()
    if (offdiagonal) call paracrystalline_SS()
!****************************************************************************************************
    

    
    !configuration counter
    complete = 0

 !****************************************************************************************************
                        !To understand these call functions you need to have a basic understanding of how mpi works. Please go through
                        !mpi tutorials available on the internet. 
 
   call MPI_INIT(ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, iproc, ierr)
   call para_range(1,conf_max ,nproc,iproc,ista,iend)
   
  !****************************************************************************************************
     
     do config = ista, iend
    
        !initialize the hamiltonian and allocate (exciton_common_local)
        call allocate_hev()
        
  !*******************************************************************************************************************************!
                        !We build the hamiltonian here; the one particle hamiltonian, two particle hamiltonian and the oneparticle -  two particle hamiltonian
       
        if ( one_state ) call hamiltonian_1p(config)
        if ( two_state ) call hamiltonian_2p(config)
        if ( one_state .and. two_state ) call hamiltonian_1p2p(config)
 !*******************************************************************************************************************************!        
        
      ! if ( LL) call long_long()
       if ( SS ) call short_short(config)
      ! if ( LS ) call long_short()
      
 !*******************************************************************************************************************************!
                        !We diagonalize the hamiltonian and find the eigen values and the eigen vectors
  
       call cpu_time(start)
       call dsyevr_diagonalize( h, kount, eval, kount2, rrange, iu ) 
       call cpu_time(finish)
        print*, 'Diagonalization Time::', finish-start
 !*******************************************************************************************************************************!        
 

 !*******************************************************************************************************************************!
                        !Here we calculate the Charge Modulation Spectra
     
        call cms_osc()
        call cms_osc_sum( sumoscx, sumoscy )
        call cms_spec( ab_x, ab_y )
 !*******************************************************************************************************************************!                
      

 !*******************************************************************************************************************************!
                        !Here we calculate the polaron coherence function. Before we calculate the coherence function we need to initialize it

    if (.not. allocated( cohfxn ) ) then 
            allocate( cohfxn( -lattice_x+1:lattice_x-1, &
                              -lattice_y+1:lattice_y-1, &
                              -lattice_z+1:lattice_z-1 ) )
            cohfxn = 0.d0
        end if        
        if (.not. allocated( mcohfxn ) ) then 
            allocate( mcohfxn( -lattice_x+1:lattice_x-1, &
                              -lattice_y+1:lattice_y-1, &
                              -lattice_z+1:lattice_z-1 ) )
            mcohfxn = 0.d0
        end if 
        !allocate the temp function and initialize to zero every time!
        if (.not. allocated( cohtmp ) )                 &
            allocate( cohtmp( -lattice_x+1:lattice_x-1, &
                              -lattice_y+1:lattice_y-1, &
                              -lattice_z+1:lattice_z-1 ) )
        cohtmp = 0.d0
       
        call cms_cohfxn( cohtmp )
        
        cohfxn = cohfxn + cohtmp
         complete = complete + 1
 !*******************************************************************************************************************************!                            
    

 end do

 !*******************************************************************************************************************************!
                        !These are mpi reduce functions. Please go though how mpireduce works

     
     Call MPI_REDUCE(bin1, mbin1, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     Call MPI_REDUCE(bin2, mbin2, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     Call MPI_REDUCE(bin3, mbin3, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     Call MPI_REDUCE(bin4, mbin4, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr) 
     Call MPI_REDUCE(bin5, mbin5, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     Call MPI_REDUCE(bin6, mbin6, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     Call MPI_REDUCE(bin7, mbin7, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     Call MPI_REDUCE(bin8, mbin8, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     Call MPI_REDUCE(bin9, mbin9, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr) 

     Call MPI_REDUCE(ab_x, mab_x, spec_step, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     ab_x = mab_x/( mbin7 * 1.d0)
     
     
     Call MPI_REDUCE(ab_y, mab_y, spec_step, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     ab_y = mab_y/( mbin7 * 1.d0)
     
     Call MPI_REDUCE(cohfxn, mcohfxn, 64, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     cohfxn = mcohfxn/( mbin7 * 1.d0 )
     !if (iproc == 0) call MPI_FINALIZE(ierr)
     
     Call MPI_REDUCE(sumoscx, msumoscx, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     sumoscx = msumoscx/( mbin7 * 1.d0 )*hw
          
     Call MPI_REDUCE(sumoscy, msumoscy, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     sumoscy = msumoscy/( mbin7 * 1.d0 )*hw
     
     if (iproc == 0) then    
     
     print*, 'mbin1', mbin1
     print*, 'mbin2', mbin2
     print*, 'mbin3', mbin3
     print*, 'mbin4', mbin4
     print*, 'mbin5', mbin5
     print*, 'mbin6', mbin6
     print*, 'mbin7', mbin7
     print*, 'mbin8', mbin8
     print*, 'mbin9', mbin9
     print*, '***************************'
     print*, ' Finishing Up...'
     print*, '***************************'

    !Normalize
    print*, 'bin9', bin9
    !ab_x = ab_x/( conf_max * 1.d0)
    !ab_y = ab_y/( conf_max * 1.d0)
    !cohfxn = cohfxn/( conf_max * 1.d0 )
    !sumoscx = sumoscx/( conf_max * 1.d0 )*hw
    !sumoscy = sumoscy/( conf_max * 1.d0 )*hw    !multiply by hw to keep in units of (cm-1)
    print*,            'Absorption Sums:'
    print*,'             x  : ', sumoscx
    print*,'             y  : ', sumoscy
    print*,'             sum: ', sumoscx + sumoscy
    write( sstring, * ) 'Absorption Sum (xysum),', sumoscx, ',', sumoscy, ',', sumoscx + sumoscy
    call para_out_append( sstring )

    !write the spectra (local)  
    call cms_out( ab_x, ab_y )
    call cohfxn_out( cohfxn )
    
    print*, ' Done'
    
end if  
    call MPI_FINALIZE(ierr)
    
end program


!****************************************************************************************************
                                                                                               !main program ends
!****************************************************************************************************