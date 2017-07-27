!***************************************************************************************
                       !Module containing all the Variables for the program
!*************************************************************************************** 

module common_variables 
implicit none
 integer conf_max, conf_process
        real*8 dwidth, ycorrlen, xcorrlen, time_start, time_finish,    &
               delta_da,sigma
        real*8, allocatable :: offset(:,:)
          !lattice numbers
    integer :: lattice_x = 1
    integer :: lattice_y = 1
    integer :: lattice_z = 1
    
    !mpi
    integer :: dim2 = 10
  
    !vibration
    real*8 :: hw = 1400.d0
    real*8 :: s  = 1.d0
    real*8 :: s_cat = 0.d0
    real*8 :: s_ani = 0.d0
    integer:: vibmax = 0

    
    !coupling
    real*8 :: jo_x = 0.d0
    real*8 :: jo_y = 0.d0
    real*8 :: jo_z = 0.d0
    real*8 :: ct_u = 0.d0
    real*8 :: ct_v = 0.d0
    real*8 :: te_x = 0.d0
    real*8 :: te_y = 0.d0
    real*8 :: te_z = 0.d0
    real*8 :: th_x = 0.d0
    real*8 :: th_y = 0.d0
    real*8 :: th_z = 0.d0

    !coupling
    logical :: extended_cpl = .false.
    real*8, allocatable :: jo_ex(:,:,:)
    character*256 extended_cpl_file

    !task title
    character*256 task_title

    !hamiltonian counters
    integer :: kount = 0
    integer :: kount_1p = 0
    integer :: kount_2p = 0
    integer :: kount_3p = 0
    integer :: kount_ct = 0
    integer :: kount_ct2 = 0
    integer :: kount_lattice = 0
    integer :: nnz = 0
    

    !indexes
    integer, allocatable :: nx_lattice(:,:,:)
    integer, allocatable :: nx_1p(:,:)
    integer, allocatable :: indx(:)
    integer, allocatable :: jndx(:)
    integer, allocatable :: nx_2p(:,:,:,:)
    real *8, allocatable :: disorder(:,:,:)  
    real, allocatable :: pdisorder(:,:,:,:,:)
    
    !the hamiltonian
    real*8, allocatable :: h(:,:)
    real*8, allocatable :: vl(:)
    real*8, allocatable :: eval(:)
    
    !franck condon factors
    real*8, allocatable :: fc_gf(:,:)
    real*8, allocatable :: fc_gc(:,:)
    real*8, allocatable :: fc_ga(:,:)
    real*8, allocatable :: fc_af(:,:)
    real*8, allocatable :: fc_cf(:,:)

    !constants
    real*8, parameter :: pi = 4.d0*datan(1.d0)
    real*8, parameter :: ev = 8065.d0
    real*8, parameter :: hc = 1.23984193d3 * ev !(plancks constant times the speed of light in nm*hw )
    real*8, parameter :: kb = 0.6956925d0        !units of cm-1 k
    real*8, parameter :: beta = 2.35
    !empty
    integer, parameter :: empty = -1

    !multiparticle states
    logical :: one_state    =.true.
    logical :: two_state    =.false.
    logical :: LL    =.true.
    logical :: SS    =.false.
    logical :: LS = .false.
    logical :: offdiagonal = .true.
    logical :: diagonal = .false.
    logical :: three_state  =.false.
    logical :: ct_state     =.false.
    logical :: ct2_state    =.false.
    
    !periodic
    logical periodic
    
    !oscillator strength
    real*8, allocatable :: osc_x(:), osc_y(:), osc_z(:)
    real*8, allocatable :: plosc_x(:,:), plosc_y(:,:), plosc_z(:,:)
    real*8, allocatable :: ux(:), uy(:), uz(:)
    real*8 :: abs_lw = 100.d0
    real*8 :: mon_tr_e = 0.d0
    integer, parameter  :: spec_step = 2600        !for 10 cm-1 resolution

    !number of threads
    integer :: numthreads = 1
    
    !dielectric
    real*8 :: dielectric=1.d0
    
    !ct states trunctation
    logical :: ct_truncate = .false.
    logical :: cty_truncate = .false.
    integer :: ct_range = 1000
    integer :: cty_range = 1000
    logical :: two_truncate = .false.
    integer :: two_range = 1000
    
    logical :: lorentzian = .false.
    logical :: linewidthprogress = .false.
    real*8 :: delta_sigma = 0.d0
    real*8, allocatable :: vibsum_evec(:), vibsum_state(:)
    logical :: abs_freq_dep = .false.
    
    
    integer :: kount2
    integer :: iu=1
    character(1) :: rrange='A'

!***************************************************************************************
          !Petsc and Slepc variables and headers required for the program
!*************************************************************************************** 

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>        
#include <petsc/finclude/petscmat.h>
#include <slepc/finclude/slepcsys.h>
#include <slepc/finclude/slepceps.h>

!*************************************************************************************** 

 
       Mat            A                          ! Matrix A
       EPS            eps                        ! Eigen Problem Solver
       EPSType        tname                      ! Type of Eigen Problem Solver  
       PetscReal      error                       
       PetscRandom    rnd                        ! Petsc Random no. generator 
       PetscInt       i, Istart, Iend, j
       Vec            xr_g,xi_g
       Vec            xr_e,xi_e
       Vec            xr_g_0, xr_e_0
       PetscScalar    kr_g,ki_g
       PetscScalar    kr_e,ki_e,value
       PetscInt       nev, maxit, its, nconv,n,mpd
       PetscMPIInt    rank,color,row_rank,row_size,new_comm
       PetscErrorCode ierr
       PetscBool      flg,terse
       PetscOffset    xx_ig,xx_ie
       double precision xx_vg(1),xx_ve(1)

    end module

        
!***************************************************************************************
                        !main program starts 
!***************************************************************************************


program cms
use common_variables
implicit none
integer ista,ien,nproc
integer config,complete


  
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
    
  
!***************************************************************************************
                       !set the parameters (local)
!***************************************************************************************    
    call set_para()
    
   
!***************************************************************************************
                       !Indexing the lattice
!***************************************************************************************    
   

    call index_lattice()
    if ( one_state )   call index_1p()
    if ( two_state )   call index_2p()    
    
    !output the kount
    print*, '***************************'
    print*, 'kount   :', kount                      !total basis set size
    print*, 'kount 1p:', kount_1p
    print*, 'kount 2p:', kount_2p
    print*, '***************************'
    
    
!***************************************************************************************
                       !Programming with Petsc
!***************************************************************************************       

  call SlepcInitialize(PETSC_NULL_CHARACTER, ierr)      ! Initializes the slepc library, SlepcInitialize() calls the PetscInitialize(), should be called in the beginning of the program
  call MPI_COMM_SIZE(PETSC_COMM_WORLD, nproc, ierr)       ! calls the no. of processors 
  call MPI_Comm_Rank(PETSC_COMM_WORLD, rank, ierr)      ! calls the rank of the processor
  color = rank/1
  call MPI_Comm_Split(PETSC_COMM_WORLD,color,rank,new_comm,ierr)
  call MPI_Comm_rank(new_comm,row_rank,ierr)
  call MPI_Comm_size(new_comm,row_size,ierr)
 
                      
   
      call MatCreate(new_comm,A,ierr)                                ! creates the matrix
      call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,kount,kount,ierr)         ! Sets the sizes
      call MatSetFromOptions(A,ierr)
      call MatSetUp(A,ierr)                                                  !Sets up the internal matrix data structures for the later use. 
     
           
!***************************************************************************************
                       !Build the hamiltonian
        
        if ( one_state ) call hamiltonian_1p_sparse(config)
        if ( two_state ) call hamiltonian_2p_sparse(config)
        if ( one_state .and. two_state ) call hamiltonian_1p2p_sparse(config)
        
        
!***************************************************************************************   
     ! Should be called after setting the values. Assembles all the matrix values.
  
      call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY,ierr)                    !  Begins assembling the matrix. This routine should be called after completing all calls to MatSetValues(). 
      call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)                       !  Completes assembling the matrix. This routine should be called after MatAssemblyBegin(). 
       
   
      

!***************************************************************************************
                       !Programming with SLEPC
!***************************************************************************************  

!   ** Create eigensolver context
      call EPSCreate(new_comm,eps,ierr)                  ! create the eigen problem solver

 !    ** Set operators. In this case, it is a standard eigenvalue problem 
      call EPSSetOperators(eps,A,PETSC_NULL_OBJECT,ierr)          
      call EPSSetProblemType(eps,EPS_HEP,ierr)                   ! My matrix is hermitian hence EPS_HEP
      nev = 1000                                               ! it is the total no. of eigen values   
      nconv = 2*nev                                              ! No. of converged eigen values
      mpd = nconv                                                
      call EPSSetDimensions(eps,nev,nconv,mpd,ierr)              ! We can set how many eigen values we need, if this command is not there, the program by default finds only one eigen pair           

 !     ** Set solver parameters at runtime
      call EPSSetFromOptions(eps,ierr)

 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 !     Solve the eigensystem
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
       call EPSSetWhichEigenpairs(eps,eps_smallest_real,ierr)     ! Sort the eigen values according to priority. If not done it store the eigen values randomly, I need from smallest to highest
     !  call EPSSetType(eps,EPSJD,ierr)
       call cpu_time(start)
       call EPSSolve(eps,ierr)                                    ! Solve
       call cpu_time(finish)

       print*, 'Diagonalization time:', finish-start 
       call EPSGetConverged (eps,nconv,ierr)                      ! Find the no. of converged eigenpairs
      if (rank .eq. 0) then
       write(*,150) nconv
     endif 
      150  format (' Number of converged eigenpairs:',I2) 
         
       

!    ** Optional: Get some information from the solver and display it
       call EPSGetType(eps,tname,ierr)
      if (rank .eq. 0) then
       write(*,120) tname
     endif
  120  format (' Solution method: ',A)
       call EPSGetDimensions(eps,nev,PETSC_NULL_INTEGER,                 &
     &                      PETSC_NULL_INTEGER,ierr)
      if (rank .eq. 0) then
      write(*,130) nev
       endif
  130  format (' Number of requested eigenvalues:',I4)
   
call EPSDestroy(eps,ierr)
call MatDestroy(A,ierr)

  
     call MPI_Comm_free(new_comm,ierr)
     call SlepcFinalize(ierr)
   
    
end program

subroutine set_para()
     use common_variables
    implicit none

    logical         exists
    character*100   buffer, label, fname
    integer         fno, ios, line, pos, errstat
    parameter       ( fno = 90 )

    !Set the default parameters
    print*, 'Setting the default parameters.'
    lattice_x       = 10
    lattice_y       = 1
    vibmax          = 0
    hw              = 1400.d0
    s               = 1.d0
    jo_x            = -0.40d0     !These are really te or th, but using exciton subroutines
    jo_y            = -0.15d0     !I have to use the jo variables
    delta_da        = 0.d0        !donor - acceptor energy difference. Set to 0 for homopolymer
    task_title      = 'test'
    abs_lw          = 350.d0
    lorentzian      = .false.
    rrange          = 'A'         !By default, find all of the eigenstates
    iu              = 1           
   ! rrange          = 'I'
    !iu              = 10
    conf_max        = 1      !disorder conf parameters
    dwidth          = 0.d0
    sigma           = 0.30d0
    xcorrlen        = 0.d0
    ycorrlen        = 0.d0
    one_state       = .true.
    two_state       = .true.
    LL = .true.
    LS = .false.
    SS = .false.
    diagonal = .false.
    offdiagonal = .false.
    two_truncate    = .false.
    two_range       = 3

      
    goto 1010
    !Read the name of the input file
    call get_command_argument(1, fname, status=errstat )
    if ( errstat .ne. 0 ) goto 1010
    inquire( file = trim(fname), exist = exists)
    if ( .not. exists ) then
        print*, 'Input file not found...aborting'
        stop
    end if
    
    !Open the file to read in
    open( unit = fno, file = fname, status = 'old', action = 'read' )
    !-----------------------------------!
    !   This reads the control file
    !   The ios changes if end of record
    !   or end of file is reached
    !-----------------------------------!
    ios = 0
    line = 0
    print*, 'Reading the input file...'
    do while ( ios == 0 )
        read( fno, '(a)', iostat=ios ) buffer
        if( ios == 0 ) then
            line = line + 1

            !Find first instance of whitespace
            pos = scan( buffer, ' ' )
            label = buffer( 1:pos )
            buffer = buffer( pos + 1:)

            select case ( label )
            case('lattice_x')
                read( buffer, *, iostat=ios ) lattice_x
                print*, '    Setting lattice_x to: ', lattice_x
            case('lattice_y')
                read( buffer, *, iostat=ios ) lattice_y
                print*, '    Setting lattice_y to: ', lattice_y
            case('vibmax')
                read( buffer, *, iostat=ios ) vibmax
                print*, '    Setting vibmax to: ', vibmax
            case('hw')
                read( buffer, *, iostat=ios ) hw
                print*, '    Setting vibrational energy to (cm-1): ', hw
            case('s')
                read( buffer, *, iostat=ios ) s
                print*, '    Setting s to : ', s
            case('te_x')
                read( buffer, *, iostat=ios) jo_x
                print*, '    Setting te_x to (eV): ', jo_x
            case('th_x')
                read( buffer, *, iostat=ios) jo_x
                print*, '    Setting th_x to (eV): ', jo_x    
            case('delta_da')
                read( buffer, *, iostat=ios) delta_da
                print*, '    Setting delta_da to (cm-1): ', delta_da           
            case('te_y')
                read( buffer, *, iostat=ios) jo_y
                print*, '    Setting te_y to (eV): ', jo_y
            case('th_y')
                read( buffer, *, iostat=ios) jo_y
                print*, '    Setting th_y to (eV): ', jo_y
            case('task_title')
                read( buffer, *, iostat=ios) task_title
                print*, '    Setting task_title to: ', trim(task_title)
            case('abs_lw')
                read( buffer, *, iostat=ios) abs_lw
                print*, '    Setting the linewidth to (cm-1): ', abs_lw
            case('two_range')
                read( buffer, *, iostat=ios) two_range
                print*, '    Setting the vibrational radius to: ', two_range
            case('two_truncate')
                read( buffer, *, iostat=ios) two_truncate
                if ( two_truncate ) print*, '    Vibrational radius is on'
                if ( .not. two_truncate ) print*, '    Vibrational radius is off'
            case('lorentzian')
                read( buffer, *, iostat=ios) lorentzian
                if ( lorentzian ) print*, '    Setting the lineshape to:', &
                                          ' lorentzian'
                if ( .not.lorentzian ) print*, '    Setting the lineshape to:', &
                                          ' gaussian'
            case('one_state')
                read( buffer, *, iostat=ios) one_state
                if ( one_state ) print*, '    One particle states are turned on.'
                if ( .not.one_state ) print*, '    One particle states are turned off.'
            case('two_state')
                read( buffer, *, iostat=ios) two_state
                if ( two_state ) print*, '    Two particle states are turned on.'
                if ( .not.two_state ) print*, '    Two particle states are turned off.'
            case('rrange')
                read( buffer, *, iostat=ios) rrange
                if ( rrange .eq. 'A' ) print*, '    dsyevr will find all  eigenvectors'
                if ( rrange .eq. 'I' ) print*, '    dsyevr will find a range of eigenvectors'
            case('iu')
                read( buffer, *, iostat=ios) iu
                print*, '    desyevr will look for the lowest', iu, ' eigenvectors'           
            case('conf_max')
                read( buffer, *, iostat=ios) conf_max
                print*, '    Setting conf_max to: ', conf_max
            case('dwidth')
                read( buffer, *, iostat=ios) dwidth
                print*, '    Setting disorder width to (cm-1): ', dwidth
            case('xcorrlen')
                read( buffer, *, iostat=ios) xcorrlen
                print*, '    Setting xcorrlen length to: ', xcorrlen            
            case('ycorrlen')
                read( buffer, *, iostat=ios) ycorrlen
                print*, '    Setting ycorrlen length to: ', ycorrlen
            case default
                print*, '    invalid label at line, ', line
            end select
        end if

    end do

    close(fno)

1010 continue
    print*, 'Calculating derived parameters.'
    !Normalize parameters to hamiltonian units of vib quanta
    jo_x = jo_x * eV / hw
    jo_y = jo_y * eV / hw
    abs_lw = abs_lw / hw
    dwidth = dwidth / hw
    delta_da = delta_da / hw
    if ( vibmax == 0 ) s = 0.d0
    
!    call omp_set_num_threads(1)

end subroutine



!**********************************************!
!                        index the lattice                                           !
!**********************************************!
subroutine index_lattice()
    use common_variables
    implicit none 

    integer lx, ly, lz
    
    if ( .not. allocated( nx_lattice ) ) then 
        allocate( nx_lattice( lattice_x, lattice_y, lattice_z ) )
    end if
    nx_lattice = empty
    
    do lx = 1, lattice_x
    do ly = 1, lattice_y
    do lz = 1, lattice_z
        kount_lattice = kount_lattice + 1
        nx_lattice( lx, ly, lz ) = kount_lattice
    end do
    end do
    end do
        
end subroutine

!**********************************************!
!                        index the 1 p states                                   !
!**********************************************!
subroutine index_1p()
    use common_variables
    implicit none 
    
    integer lx, ly, lz, lxyz, vib

    !allocate the index
    if ( .not. allocated( nx_1p ) ) then 
            allocate( nx_1p ( kount_lattice, 0:vibmax ) )
    end if
    nx_1p = empty
    
    do lx = 1, lattice_x
    do ly = 1, lattice_y
    do lz = 1, lattice_z
    do vib = 0, vibmax
        kount = kount + 1
        kount_1p = kount_1p + 1
        lxyz = nx_lattice( lx, ly, lz )
        nx_1p( lxyz, vib ) = kount
        write (9,*) 'nx_1p(', lxyz, vib,' )', kount
!                print*, lxyz, vib, kount
    end do
    end do
    end do
    end do
        
end subroutine
!**********************************************!
!                        index the 2 p states                                   !
!**********************************************!

subroutine index_2p()
    use common_variables
    implicit none 
    
    integer lx, ly, lz, lxyz, vib
    integer lxv, lyv, lzv, lxyzv, vibv

    !allocate the index
    if ( .not. allocated( nx_2p ) ) then 
            allocate( nx_2p ( kount_lattice, 0:vibmax, kount_lattice, 1:vibmax ) )
    end if
    nx_2p = empty
    
    do lx = 1, lattice_x
    do ly = 1, lattice_y
    do lz = 1, lattice_z
    do vib = 0, vibmax
        do lxv = 1, lattice_x
        do lyv = 1, lattice_y
        do lzv = 1, lattice_z
!                do vib = 0, vibmax
        do vibv = 1, vibmax
            lxyz = nx_lattice( lx, ly, lz )
            lxyzv = nx_lattice( lxv, lyv, lzv ) 
            if ( lxyz == lxyzv ) cycle                !not a 2 p state
            if (vibv + vib > vibmax) cycle
            if ( two_truncate .and. two_range*1.d0 < dsqrt( 1.d0*((lx-lxv)**2+(ly-lyv)**2 + (lzv-lzv)**2)) ) cycle      !truncate at two_range
            kount = kount + 1
            kount_2p = kount_2p + 1
            nx_2p( lxyz, vib, lxyzv, vibv ) = kount
            write (9,*) 'nx_2p(', lxyz, vib, lxyzv, vibv,' )', kount
            !write (7,*) 'nx_2p(',lxyz,vib,lxyzv,vibv,')', kount
!                        print*, lxyz, vib,lxyzv, vibv, kount
        end do
        end do
        end do
        end do
    end do
    end do
    end do
    end do
        
end subroutine

!**********************************************************************************!
!                        create 1p sparse hamiltonian                                 !
!**********************************************************************************!
subroutine hamiltonian_1p_sparse(config)
    use common_variables
    implicit none
    
    integer lx1, ly1, lz1, lxyz1, vib1, hx1
    integer lx2, ly2, lz2, lxyz2, vib2, hx2,config
    real*8 get_j,coupling,val,dis_order
    real*8 SS_disorder
    if ( .not. allocated( vibsum_state ) ) allocate( vibsum_state(kount) )
    vibsum_state = 0.d0
 call MatGetOwnershipRange(A,Istart,Iend,ierr)
    do lx1 = 1, lattice_x
    do ly1 = 1, lattice_y
    do lz1 = 1, lattice_z
    do vib1 = 0, vibmax
    lxyz1 = nx_lattice( lx1, ly1, lz1 )
      
     
      
       hx1 = nx_1p( lxyz1, vib1 )
       if ( hx1 == empty ) cycle

       hx1 = hx1-1
       if (hx1 >= Istart .and. hx1<Iend) then
       
        !diagonal
        val = vib1*1.d0
      !  if (diagonal) then
 !print*, 'vib1', vib1
 !print*, 'disorder (', config,lx1,ly1,')', disorder(config,lx1,ly1)
  !      SS_disorder = disorder(config,lx1,ly1)
 !print*, 'SS_disorder', SS_disorder
   !     dis_order = val + SS_disorder
!end if
      if (val .eq. 0) then
continue
else
 call MatSetValues ( A, 1, hx1, 1, hx1, val,INSERT_VALUES,ierr)        !Set the values
  end if
     
        do lx2 = 1, lattice_x
        do ly2 = 1, lattice_y
        do lz2 = 1, lattice_z
        do vib2 = 0, vibmax
        lxyz2 = nx_lattice( lx2, ly2, lz2 )
        hx2 = nx_1p( lxyz2, vib2 )
            if ( hx2 == empty ) cycle
            hx2=hx2-1
            if ( hx2 == hx1 ) cycle
     
  
   
            !off diagonal
     coupling = 1
 if (coupling .eq. 0)  then
continue
else
     call MatSetValues ( A, 1,hx1, 1, hx2, coupling ,INSERT_VALUES,ierr)
    end if         
        


       
        end do
        end do
        end do
        end do
end if
  
    end do
    end do
    end do
    end do

end subroutine
  
!**********************************************************************************!
!                        create 1p2p sparse hamiltonian                            !
!**********************************************************************************!

subroutine hamiltonian_1p2p_sparse(config)
    use common_variables
    implicit none
    
    integer lx1, ly1, lz1, lxyz1, vib1, hx1,config
    integer lx2, ly2, lz2, lxyz2, vib2, hx2
    integer lx2v, ly2v, lz2v, lxyz2v, vib2v
    real*8 get_j,val,coupling
    
    call MatGetOwnershipRange(A,Istart,Iend,ierr)
    
    do lx1 = 1, lattice_x
    do ly1 = 1, lattice_y
    do lz1 = 1, lattice_z
    do vib1 = 0, vibmax
        lxyz1 = nx_lattice( lx1, ly1, lz1 )
        hx1 = nx_1p( lxyz1, vib1 )
        if ( hx1 == empty ) cycle
        hx1 = hx1-1
        
        if (hx1 >= Istart .and. hx1<Iend) then  
      
        do lx2 = 1, lattice_x
        do ly2 = 1, lattice_y
        do lz2 = 1, lattice_z
        do vib2 = 0, vibmax
            lxyz2 = nx_lattice( lx2, ly2, lz2 )
            
            !these are the only ones that are non-zero
            lx2v = lx1 
            ly2v = ly1
            lz2v = lz1
            lxyz2v = nx_lattice( lx2v, ly2v, lz2v )
            do vib2v = 1, vibmax

                hx2 = nx_2p( lxyz2, vib2, lxyz2v, vib2v )
                if ( hx2 == empty ) cycle
                hx2=hx2-1  
                
              coupling = 1
           if (coupling .eq. 0)  then
             continue
             else
     call MatSetValues ( A, 1,hx1, 1, hx2, coupling ,INSERT_VALUES,ierr)
     call MatSetValues ( A, 1,hx2, 1, hx1, coupling ,INSERT_VALUES,ierr)
            end if   
               
               
           
            end do
        end do
        end do
        end do
        end do
    
    end if

    end do
    end do
    end do
    end do
   
end subroutine


!**********************************************************************************!
!                        create 2p sparse hamiltonian                              !
!**********************************************************************************!


subroutine hamiltonian_2p_sparse(config)
    use common_variables
    implicit none
    
    integer lx1, ly1, lz1, lxyz1, vib1, hx1,config
    integer lx1v, ly1v, lz1v, lxyz1v, vib1v
    integer lx2, ly2, lz2, lxyz2, vib2, hx2
    integer lxyz2v, vib2v
    real*8 get_j,coupling,val
    real*8 SS_disorder,dis_order
    call MatGetOwnershipRange(A,Istart,Iend,ierr) 

    do lx1 = 1, lattice_x
    do ly1 = 1, lattice_y
    do lz1 = 1, lattice_z
    do vib1 = 0, vibmax
        lxyz1 = nx_lattice( lx1, ly1, lz1 )
    
        do lx1v = 1, lattice_x
        do ly1v = 1, lattice_y
        do lz1v = 1, lattice_z
        do vib1v = 1, vibmax
            lxyz1v = nx_lattice( lx1v, ly1v, lz1v )
    
            
            hx1 = nx_2p( lxyz1, vib1, lxyz1v, vib1v )
            
           
            if ( hx1 == empty ) cycle
            hx1=hx1-1
            if (hx1 >= Istart .and. hx1<Iend) then
            val = (vib1 +vib1v)*1.d0
          ! if (diagonal) then
       ! SS_disorder = disorder(config,lx1,ly1)
! print*, 'SS_disorder', SS_disorder
       ! dis_order = val + SS_disorder
!end if
            if (val .eq. 0) then
            continue
            else
            call MatSetValues ( A, 1, hx1, 1, hx1, val,INSERT_VALUES,ierr)
            end if         
         
            do lx2 = 1, lattice_x
            do ly2 = 1, lattice_y
            do lz2 = 1, lattice_z
                lxyz2 = nx_lattice( lx2, ly2, lz2 )
            do vib2 = 0, vibmax
            
                !linker type
                lxyz2v = lxyz1v
                vib2v = vib1v
                hx2 = nx_2p( lxyz2, vib2, lxyz2v, vib2v )
                hx2=hx2-1
                if ( hx2 == empty .or. hx2 == hx1 ) then
                    continue
                else
               coupling = 1
              if (coupling .eq. 0)  then
              continue
              else
              call MatSetValues ( A, 1,hx1, 1, hx2, coupling ,INSERT_VALUES,ierr)
              end if
                   
                end if
        
                !exchange type
                if ( lxyz1v == lxyz2 ) then
                    lxyz2v = lxyz1
                    do vib2v = 1, vibmax
                        hx2 = nx_2p(  lxyz2, vib2, lxyz2v, vib2v )
                        if ( hx2 == empty) cycle
                        hx2=hx2-1
                        if (hx2 == hx1 ) cycle
                        
                        coupling = 1
                  if (coupling .eq. 0)  then
                      continue
                          else
                   call MatSetValues ( A, 1,hx1, 1, hx2, coupling ,INSERT_VALUES,ierr)
                      end if   
                     
                    end do
               end if
            end do
            end do
            end do
            end do
    
        end if
        
        end do
        end do
        end do
        end do
    end do
    end do
    end do
    end do
             
end subroutine

