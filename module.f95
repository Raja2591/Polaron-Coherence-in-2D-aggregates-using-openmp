!******************************************************************************
!This module contains all the variables that the program needs
!******************************************************************************


module common_variables 
implicit none
 integer conf_max, conf_process
        real*8 dwidth, ycorrlen, xcorrlen, time_start, time_finish,    &
               delta_da
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

    !doping
    real*8 :: dintra = 0.d0
    real*8 :: dinter = 0.d0
    real*8 :: r = 0.d0
 
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

    !indexes
   integer, allocatable :: nx_lattice(:,:,:)
    integer, allocatable :: nx_1p(:,:)
    integer, allocatable :: nx_2p(:,:,:,:)
    real *8, allocatable :: disorder(:,:,:)
    real *8, allocatable :: Force(:,:,:) 
    real, allocatable :: pdisorder(:,:,:,:,:)  
    real *8, allocatable :: te (:,:,:)
    real *8, allocatable :: D(:,:,:) 
    
    !the hamiltonian
    real*8, allocatable :: h(:,:)
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
    real*8, parameter :: angstrom = 0.0000000001 !units of m
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
    logical :: dopant = .false.
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
    end module


! Module ends