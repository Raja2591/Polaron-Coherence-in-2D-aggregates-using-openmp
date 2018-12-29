!***************************************************************************************!
!    This program calculates the charge modulation spectra for a 2 dimensional
!    aggregate. It is able to run over multiple processors using OpenMP allowing
!    for larger system sizes. 
!       
!                                                   Author: Raja Ghosh
!                                                    Date  : 10/07/2018
!***************************************************************************************!

!***************************************************************************************!
!    Module contains common variables to this program
!***************************************************************************************!
module cms_common_disorder
        implicit none

        integer conf_max
        real*8 dwidth, ycorrlen, xcorrlen, time_start, time_finish,    &
               delta_da
        real*8, allocatable :: offset(:,:)
end module
!***************************************************************************************!
!    Main program
!***************************************************************************************!
program cms_polaron()
    use exciton_common_local
    use cms_common_disorder
    implicit none
    
    integer config, complete
    real*8  start, finish
    real*8 :: ab_x(spec_step) = 0.D0
    real*8 :: ab_y(spec_step) = 0.D0
    real*8 :: sumoscx = 0.d0, sumoscy = 0.d0
    real*8, allocatable :: cohfxn(:,:,:), cohtmp(:,:,:)
    integer, external   :: omp_get_thread_num
    character*128 sstring
    
    !set the parameters (local)
    call set_para()
    
    !enter simulation directory
    call mkdir( trim(task_title) )
    
    !do the indexing (exciton_common_local)
    call index_lattice()
    if ( one_state )   call index_1p()
    if ( two_state )   call index_2p()    
    
    !output the kount
    print*, '***************************'
    print*, 'kount   :', kount
    print*, 'kount 1p:', kount_1p
    print*, 'kount 2p:', kount_2p
    print*, '***************************'
    
    !make the fc table (exciton_common_local)
    call set_fctable()
        
    !set the disorder table (local)
    call set_cms_disorder_table()

    !output the parameters (local)
    call para_out_cms()
    
    !configuration counter
    complete = 0
    !$OMP PARALLEL DO PRIVATE( config, cohtmp ) SHARED( cohfxn ) REDUCTION( +:ab_x, ab_y, sumoscx, sumoscy )
    do config = 1, conf_max
    
        !initialize the hamiltonian and allocate (exciton_common_local)
        call allocate_hev()
        
        !build the hamiltonian (exciton_common_local)
        if ( one_state ) call hamiltonian_1p()
        if ( two_state ) call hamiltonian_2p()
        if ( one_state .and. two_state ) call hamiltonian_1p2p()
        call hamiltonian_cms_d(config)     !(local)
        call hamiltonian_cms_da()          !add the energy difference for the hetero polymer

        !find the eigenspectrum
        call cpu_time(start)
        call dsyevr_diagonalize( h, kount, eval, kount2, rrange, iu )    !(exciton_common_local)
        call cpu_time(finish)
!        print*, 'Diagonalization Time::', finish-start
        
        !calculate the spectra (local)
        call cms_osc()
        call cms_osc_sum( sumoscx, sumoscy )
        call cms_spec( ab_x, ab_y )
        
        !calculate the coherence function    
        !allocate the coherence function and initialize it.
        if (.not. allocated( cohfxn ) ) then 
            allocate( cohfxn( -lattice_x+1:lattice_x-1, &
                              -lattice_y+1:lattice_y-1, &
                              -lattice_z+1:lattice_z-1 ) )
            cohfxn = 0.d0
        end if        
        !allocate the temp function and initialize to zero every time!
        if (.not. allocated( cohtmp ) )                 &
            allocate( cohtmp( -lattice_x+1:lattice_x-1, &
                              -lattice_y+1:lattice_y-1, &
                              -lattice_z+1:lattice_z-1 ) )
        cohtmp = 0.d0
       
        call cms_cohfxn( cohtmp )
        
        !Only allow one thread to access cohfxn at a time
        !$OMP ATOMIC
        cohfxn = cohfxn + cohtmp
        
        complete = complete + 1
        if ( mod( complete, omp_get_num_threads() ) == 0 ) then
            print*, 'completed : ', complete, '/', conf_max
        end if
        
    end do
    !$OMP END PARALLEL DO
        
    print*, '***************************'
    print*, ' Finishing Up...'
    print*, '***************************'

    !Normalize
    ab_x = ab_x/( conf_max * 1.d0)
    ab_y = ab_y/( conf_max * 1.d0)
    cohfxn = cohfxn/( conf_max * 1.d0 )
    sumoscx = sumoscx/( conf_max * 1.d0 )*hw
    sumoscy = sumoscy/( conf_max * 1.d0 )*hw    !multiply by hw to keep in units of (cm-1)
    print*,            'Absorption Sums:'
    print'(a, f14.4)', '             x  : ', sumoscx
    print'(a, f14.4)', '             y  : ', sumoscy
    print'(a, f14.4)', '             sum: ', sumoscx + sumoscy
    write( sstring, * ) 'Absorption Sum (xysum),', sumoscx, ',', sumoscy, ',', sumoscx + sumoscy
    call para_out_append( sstring )

    !write the spectra (local)  
    call cms_out( ab_x, ab_y )
    call cohfxn_out( cohfxn )
    
    print*, ' Done'

end program
!***************************************************************************************!
!    Set the parameters
!***************************************************************************************!
subroutine set_para()
    use exciton_common_local
    use cms_common_disorder
    implicit none

    logical         exists
    character*100   buffer, label, fname
    integer         fno, ios, line, pos, errstat, num_threads
    parameter       ( fno = 90 )

    !Set the default parameters
    print*, 'Setting the default parameters.'
    lattice_x       = 4
    lattice_y       = 4
    vibmax          = 4
    hw              = 1400.d0
    s               = 1.d0
    jo_x            = -0.3d0        !These are really te or th, but using exciton subroutines
    jo_y            = -0.15d0       !I have to use the jo variables
    delta_da        = 2500.d0        !donor - acceptor energy difference. Set to 0 for homopolymer
    task_title      = 'test'
    abs_lw          = 350.d0
    lorentzian      = .false.
    rrange          = 'A'           !By default, find all of the eigenstates
    iu              = 1             !disorder conf parameters
    conf_max        = 100
    dwidth          = 2500.d0
    xcorrlen        = 0.d0
    ycorrlen        = 0.d0
    one_state       = .true.
    two_state       = .true.
    two_truncate    = .false.
    two_range       = 2

    call omp_set_num_threads( 10 )

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
            case('num_threads')
                read( buffer, *, iostat=ios) num_threads
                print*, '    Setting number of threads to: ', num_threads
                call omp_set_num_threads( num_threads )
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

end subroutine
!***************************************************************************************!
!    Write the parameters to a file
!***************************************************************************************!
subroutine para_out_cms()
    use exciton_common_local
    use cms_common_disorder
    implicit none
    
    character*64 fname
    integer fno
    parameter( fno = 16 )
    character*8 ddd
    character*10 ttt
    
    !create the file
    fname = trim(task_title)//'_para.csv'
    open( unit = fno, file = fname, action='write' )
    call date_and_time( ddd, ttt )
    
    write( fno, * ) 'parameter file for, cms_local_disorder.f95'
    write( fno, * ) 'task title, ', trim(task_title)
    write( fno, * ) 'date and time, ', ddd, ' ', ttt
    write( fno, * ) 'x dimension, ', lattice_x
    write( fno, * ) 'y dimension, ', lattice_y
    write( fno, * ) 'vibmax, ', vibmax
    write( fno, * ) 'vib energy (cm-1), ', hw 
    write( fno, * ) 'huang-rhys, ', s
    write( fno, * ) 'th x (eV), ', jo_x*hw/eV
    write( fno, * ) 'th y (eV), ', jo_y*hw/eV    
    write( fno, * ) 'abs lw (cm-1), ', abs_lw * hw
    if ( lorentzian ) write( fno, * ) 'lineshape, lorentzian'
    if ( .not. lorentzian ) write( fno, * ) 'lineshape, gaussian'    
    write( fno, * ) 'conformations, ', conf_max
    write( fno, * ) 'disorder width (cm-1), ', dwidth*hw
    write( fno, * ) 'x correlation length, ', xcorrlen
    write( fno, * ) 'y correlation length, ', ycorrlen
    write( fno, * ) 'one state, ', one_state
    write( fno, * ) 'two state, ', two_state
    write( fno, * ) 'two truncate, ', two_truncate
    write( fno, * ) 'two range, ', two_range   
    
    close( fno )
    
end subroutine
!***************************************************************************************!
!    Calculate the CMS oscillator strength
!***************************************************************************************!
subroutine cms_osc()
    use exciton_common_local
    implicit none
    
    integer ground
    parameter ( ground = 1 )
    integer lx, ly, lz, lxyz, vib, hx,  &
            lxv, lyv, lzv, lxyzv, vibv, &
            excited
    real*8  tmp_x, tmp_y
    
    If ( .not. allocated( osc_x ) ) allocate( osc_x(2:kount2), osc_y(2:kount2) )
    
    !assume the stack is linear... the transition dipole moment operator is
    !                              u=r(lx,ly)|lx,ly><lx,ly|
    !go over all excited states found during diagonalization
    do excited = 2, kount2
        
        !initialize the tmp variables
        tmp_x = 0.d0
        tmp_y = 0.d0
    
        !go over all 1p basis states
        if ( one_state ) then
        do lx = 1, lattice_x
        do ly = 1, lattice_y
        do lz = 1, lattice_z
        do vib = 0, vibmax
            lxyz = nx_lattice( lx, ly, lz )
            hx = nx_1p( lxyz, vib )
            if ( hx == empty ) cycle
            
            !set the transition dipole moments
            tmp_x = tmp_x + lx * h( hx, ground ) * h( hx, excited )
            tmp_y = tmp_y + ly * h( hx, ground ) * h( hx, excited )

        end do
        end do
        end do
        end do
        end if
        
        !go over all 2p basis states
        if ( two_state ) then
        do lx = 1, lattice_x
        do ly = 1, lattice_y
        do lz = 1, lattice_z
        do vib = 0, vibmax
            lxyz = nx_lattice( lx, ly, lz )
        do lxv = 1, lattice_x
        do lyv = 1, lattice_y
        do lzv = 1, lattice_z
        do vibv = 1, vibmax
            lxyzv = nx_lattice( lxv, lyv, lzv )
            hx = nx_2p( lxyz, vib, lxyzv, vibv )
            if ( hx == empty ) cycle
            
            !set the transition dipole moments
            tmp_x = tmp_x + lx * h( hx, ground ) * h( hx, excited )
            tmp_y = tmp_y + ly * h( hx, ground ) * h( hx, excited )
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        end if
        
        !square to get the oscillator strength
        osc_x( excited ) = tmp_x * tmp_x
        osc_y( excited ) = tmp_y * tmp_y
        
    end do
end subroutine
!***************************************************************************************!
!    Calculate the CMS oscillator strength sums
!***************************************************************************************!
subroutine cms_osc_sum( sumx, sumy )
    use exciton_common_local
    implicit none
    
    real*8, intent(inout) :: sumx, sumy
    integer excited, ground
    parameter ( ground = 1 )
    real*8 trans_e
    
    do excited = 2, kount2
        trans_e = eval( excited ) - eval( ground )
        sumx = sumx + trans_e * osc_x( excited )
        sumy = sumy + trans_e * osc_y( excited )   
    end do

end subroutine
!***************************************************************************************!
!    Calculate the CMS spectrum
!***************************************************************************************!
subroutine cms_spec( ab_x, ab_y )
    use exciton_common_local
    implicit none
    
    real*8, intent(out) :: ab_x( spec_step ), ab_y( spec_step )
    integer ground
    parameter ( ground = 1 )
    integer spec_point, excited
    real*8 spec_start, spec_end, energy, trans_e, &
           lineshape, gamma
               
    spec_start = 0.d0                              !start
    spec_end = spec_start + 10000.D0/hw            !end
      
    do spec_point = 1, spec_step 
        energy = spec_start + spec_point/( spec_step * 1.D0 )*( spec_end - spec_start )
            
        do excited = 2, kount2
            trans_e = eval( excited ) - eval( ground )
                                            
            gamma = abs_lw
            
            if ( lorentzian ) then
                lineshape = gamma/( (energy - trans_E)**2 + gamma**2 )/pi
!                lineshape = lineshape/(1.D0*kount_lattice)
            else	        !Gaussian
                lineshape = dexp( - ( ( energy - trans_E ) / gamma ) ** 2 )
                lineshape = lineshape / ( dsqrt(pi) * gamma )            !Keep gaussian normalized
!                lineshape = lineshape / ( dsqrt(pi) * gamma * kount_lattice)
            end if
 
            !cms is frequency dependent, so must multiply by the transition energy
            lineshape = lineshape*trans_e
 
            ab_x( spec_point ) = ab_x( spec_point ) + lineshape * osc_x( excited )
            ab_y( spec_point ) = ab_y( spec_point ) + lineshape * osc_y( excited )
 
        end do        
    end do
  
end subroutine
!***************************************************************************************!
!    Write the CMS spectrum
!***************************************************************************************!
subroutine cms_out( ab_x, ab_y )
    use exciton_common_local
    implicit none
    
    character(256) f_name
    integer f_no, spec_point
    real*8 spec_start, spec_end, energy,                                 &
           ab_x( spec_step ), ab_y( spec_step )
        
    !File info
    f_name = trim(task_title)//'_cms.dat'
    f_no = 2
        
    !Open the file and put headers
    open( unit = f_no, file = f_name )
    write(f_no, *) 'Task Title: ', trim(task_title)
    write(f_no, *) 'CMS Data'
    write(f_no, *) 'Energy CMS CMS '
    write(f_no, *) ' cm\+(-1) a.u. a.u. a.u. '
    write(f_no, *) ' n.a. sum x y'
        
    spec_start = 0.d0                              !start
    spec_end = spec_start + 10000.D0/hw            !end
    
    do spec_point = 1, spec_step 
        energy = spec_start + spec_point/( spec_step * 1.D0 )*( spec_end - spec_start )
        !multiply all by hw to get in units of cm-1
        write( f_no, '(5f14.7)' ) energy*hw,                                           &
                                  (ab_x( spec_point ) + ab_y( spec_point )), &
                                  ab_x( spec_point ) , ab_y( spec_point )
    end do
    close( f_no )

 End Subroutine

!***************************************************************************************!
!    Set the disorder table for use in the parallel region
!***************************************************************************************!
subroutine set_cms_disorder_table()
    use exciton_common_local
    use cms_common_disorder
    implicit none
    
    character*128 string_out
    integer vx, vy, vz, vxyz,     &
            config, idist, i, j
    parameter ( idist = 3 )                       !draw from standard normal distribution
    integer :: iseed(4) = (/47, 3093,1041,77/)    !random number seed
    real*8 a( kount_lattice, kount_lattice ), &
           aeval( kount_lattice ), &
           deltabar( kount_lattice ), &
           average, stddev
    real*8, external :: dlarnd
           
    print*, 'Setting the disorder table'
    
    !allocate ofset matrix
    if ( .not. allocated( offset ) ) allocate( offset( conf_max, kount_lattice ) )
    offset = 0.d0
    
    !keep numbers from blowing up
    if ( xcorrlen == 0.d0 ) xcorrlen = 1.d-8
    if ( ycorrlen == 0.d0 ) ycorrlen = 1.d-8
    if ( dwidth   == 0.d0 ) dwidth   = 1.d-8
    
    !correlation matrices (returns diagonalizing matrix of a inverse, and the diagonal a inverse )
    call set_cms_correlation_matrix( a, aeval )
    
    !define the disorder matrix
    do config = 1, conf_max

        !Draw the frequencies corresponding to the diagonalized correlation matrix
        !from a normal distribution. The standard deviation of each frequency is
        !given by the inverse square root of the diagonalized correlation matrix 
        do vx = 1, lattice_x    !these variables are not local!!! I should change the notation here
        do vy = 1, lattice_y
        do vz = 1, lattice_z
            vxyz = nx_lattice( vx, vy, vz )
            deltabar( vxyz ) = dlarnd( idist, iseed ) / dsqrt( aeval( vxyz ) )
        end do
        end do
        end do
        
        !now define the local disorder by back transforming the frequencies
        !from the diagonal correlation basis to the local basis
        offset( config, : ) = matmul( a, deltabar )

    end do

    
    !calculate the statistics
    average = 0.d0
    stddev = 0.d0
    do config = 1, conf_max
    do i = 1, kount_lattice
        average = average + offset( config, i )
    end do
    end do

!    do i = 1, lattice_x
!        print'(4f14.2)', ( hw*offset( 1, nx_lattice( i, j, 1 ) ), j = 1, lattice_y )
!    end do

    average = average/(conf_max*kount_lattice*1.d0)
    !now standard deviaton
    do config = 1, conf_max
    do i = 1, kount_lattice
        stddev = stddev + ( average - offset( config, i ) )**2
    end do
    end do
    stddev = dsqrt(stddev / ( conf_max*kount_lattice*1.d0 ) )
    
    print*, '********************************'
    print'(2(a8, f8.2))', 'ave: ', average*hw, ' stddev: ', stddev*hw
    print*, '********************************'
    
    write( string_out, * ) 'ave (cm-1): ,', average*hw, '  ,stddev (cm-1): ,', stddev*hw
    call para_out_append( string_out )
    
end subroutine
!***************************************************************************************!
!    Set the correlation matrix
!***************************************************************************************!
subroutine set_cms_correlation_matrix( a, aeval )
    use exciton_common_local
    use cms_common_disorder
    implicit none
    
    real*8, intent(inout) :: a(kount_lattice, kount_lattice)
    real*8, intent(out)   :: aeval( kount_lattice )
    real*8  work( kount_lattice, kount_lattice )
    integer ipiv( kount_lattice ), info, lx, ly, lz, &
            lxyz, lx2, ly2, lz2, lxyz2,              &
            dx, dy, dz
            
    !make the matrix a
    do lx = 1, lattice_x
    do ly = 1, lattice_y
    do lz = 1, lattice_z
        lxyz = nx_lattice( lx, ly, lz )
        
        do lx2 = 1, lattice_x
        do ly2 = 1, lattice_y
        do lz2 = 1, lattice_z
            lxyz2 = nx_lattice( lx2, ly2, lz2 )
            
            dx = dabs( 1.d0*(lx2-lx) ) 
            dy = dabs( 1.d0*(ly2-ly) )
            dz = dabs( 1.d0*(lz2-lz) )
            
            !is this right? i think so...it reduces to exponential distance dependence when xcorrlen==ycorrlen            
            a( lxyz, lxyz2 ) = dwidth * dwidth *                     &
                               ( dexp( -1.d0*dsqrt(                  &
                                     ( dx*dx / xcorrlen**2 +         &
                                       dy*dy / ycorrlen**2 ) ) ) )     
        end do
        end do
        end do
    end do
    end do
    end do
    
    !find the inverse to return
    !first find the LU decomposition to pass to dgetri
    call dgetrf( kount_lattice, kount_lattice, a, kount_lattice, ipiv, info )
    !now find the inverse
    call dgetri( kount_lattice, a, kount_lattice, ipiv, work, kount_lattice, info )
    
    !now find eigenvalues and eigenvectors and return
    call dev_diagonalize( a, kount_lattice, aeval )
 
end subroutine
!***************************************************************************************!
!    Build the disorder hamiltonian
!***************************************************************************************!
subroutine hamiltonian_cms_d( config )
    use exciton_common_local
    use cms_common_disorder
    implicit none
    
    integer, intent(in) :: config
    integer  lx, ly, lz, lxyz,     &
             lx2, ly2, lz2, lxyz2, &
             vib, vib2, hx
             
    !one particle states
    if ( one_state ) then
        do lx = 1, lattice_x
        do ly = 1, lattice_y
        do lz = 1, lattice_z
            lxyz = nx_lattice( lx, ly, lz )
        do vib = 0, vibmax
            hx = nx_1p( lxyz, vib )
            if ( hx == empty ) cycle
            
            !add the disorder
            h( hx, hx ) = h( hx, hx ) + offset( config, lxyz )
        end do
        end do
        end do
        end do
    end if
    
    !two particle states
    if ( two_state ) then
        do lx = 1, lattice_x
        do ly = 1, lattice_y
        do lz = 1, lattice_z
            lxyz = nx_lattice( lx, ly, lz )
        do vib = 0, vibmax
            do lx2 = 1, lattice_x
            do ly2 = 1, lattice_y
            do lz2 = 1, lattice_z
                lxyz2 = nx_lattice( lx2, ly2, lz2 )
            do vib2 = 1, vibmax
                hx = nx_2p( lxyz, vib, lxyz2, vib2 )
                if ( hx == empty ) cycle
                
                h( hx, hx ) = h( hx, hx ) + offset( config, lxyz )
            
            end do
            end do
            end do
            end do
        end do
        end do
        end do
        end do
    end if
    
end subroutine
!***************************************************************************************!
!    Build the hamiltonian for the hetero polymer
!***************************************************************************************!
subroutine hamiltonian_cms_da()
    use exciton_common_local
    use cms_common_disorder
    implicit none
    
    integer  lx, ly, lz, lxyz,     &
             lx2, ly2, lz2, lxyz2, &
             vib, vib2, hx
             
    !one particle states
    if ( one_state ) then
        do lx = 1, lattice_x
        do ly = 1, lattice_y
        do lz = 1, lattice_z
            lxyz = nx_lattice( lx, ly, lz )
        do vib = 0, vibmax
            hx = nx_1p( lxyz, vib )
            if ( hx == empty ) cycle
            
            !add the donor-acceptor offset
            if ( ( mod( lx, 2 ) == 0 .and. mod( ly, 2 ) == 0 ) .or.     &    !(even row, even column)
                 ( mod( lx, 2 ) == 1 .and. mod( ly, 2 ) == 1 ) ) then        !(odd row, odd column) this will make a checkerboard da pattern
                    h( hx, hx ) = h( hx, hx ) + delta_da
            end if
        end do
        end do
        end do
        end do
    end if
    
    !two particle states
    if ( two_state ) then
        do lx = 1, lattice_x
        do ly = 1, lattice_y
        do lz = 1, lattice_z
            lxyz = nx_lattice( lx, ly, lz )
        do vib = 0, vibmax
            do lx2 = 1, lattice_x
            do ly2 = 1, lattice_y
            do lz2 = 1, lattice_z
                lxyz2 = nx_lattice( lx2, ly2, lz2 )
            do vib2 = 1, vibmax
                hx = nx_2p( lxyz, vib, lxyz2, vib2 )
                if ( hx == empty ) cycle
                
                if ( ( mod( lx, 2 ) == 0 .and. mod( ly, 2 ) == 0 ) .or.     &    !(even row, even column)
                 ( mod( lx, 2 ) == 1 .and. mod( ly, 2 ) == 1 ) ) then            !(odd row, odd column) this will make a checkerboard da pattern
                
                    h( hx, hx ) = h( hx, hx ) + delta_da
                end if
            
            end do
            end do
            end do
            end do
        end do
        end do
        end do
        end do
    end if
    
end subroutine
!***************************************************************************************!
!    Make a directory for the run
!***************************************************************************************!
subroutine mkdir( dir )
    implicit none
        
    character(*), intent(in) :: dir
    integer fno
    character(100) fname
    character*100 inp
    parameter ( fno = 15 )
    fname=trim(dir)//'mkdir.bat'

    open( unit = fno, file = fname )
    write( fno, * ) '@echo off'
    write( fno, * ) 'if not exist '//dir//' mkdir '//dir
    call get_command_argument( 1, inp )
    write( fno, * ) 'move  /Y '//trim(inp)//' '//trim(dir)//' 1>NUL 2>&1'
    close (fno)
    
    call system( fname )
    
    open( unit = fno, file = fname )
    close( fno, status = 'delete' )
    
    print*, 'entering \'//dir
    call chdir( dir )
end subroutine
!***************************************************************************************!
!    Calculate the coherence function, here it is
!    C(r)=<Psi|sumn d_n^dagger d_(n+r)|Psi>
!    This seems to be behaving ok, but I am not sure that it is right!!!!!!!!
!    Check with Spano.
!***************************************************************************************!
subroutine cms_cohfxn( cohfxn )
    use exciton_common_local
    implicit none
    
    real*8, intent(inout) :: cohfxn( -lattice_x+1:lattice_x-1, &
                                     -lattice_y+1:lattice_y-1, &
                                     -lattice_z+1:lattice_z-1 )
    integer lx, ly, lz, lxyz,     &
            rx, ry, rz,           &
            lx2, ly2, lz2, lxyz2, &
            vib, vib2, ground,    &
            hx, hx2,              &
            lxv, lyv, lzv, lxyzv, &
            vibv, vibg,           &
            lxyz2v, vib2v
    parameter ( ground = 1, vibg = 0 )

    !choose r
    do rx = -lattice_x + 1, lattice_x - 1
    do ry = -lattice_y + 1, lattice_y - 1
    do rz = -lattice_z + 1, lattice_z - 1
    
        !one-particle one-particle
        if ( one_state ) then
        do lx = max( 1-rx, 1 ), min( lattice_x - rx, lattice_x )
        do ly = max( 1-ry, 1 ), min( lattice_y - ry, lattice_y )
        do lz = max( 1-rz, 1 ), min( lattice_z - rz, lattice_z )
            lx2 = lx + rx
            ly2 = ly + ry
            lz2 = lz + rz
            lxyz = nx_lattice( lx, ly, lz )
            lxyz2 = nx_lattice( lx2, ly2, lz2 )
            
            !if on the same molecule, vib must be the same on both
            if ( lxyz == lxyz2 ) then
                do vib = 0, vibmax
                    hx = nx_1p( lxyz, vib )
                    if ( hx == empty ) cycle
                    cohfxn( rx, ry, rz ) = cohfxn( rx, ry, rz ) +                 &
                                           h( hx, ground ) * h( hx, ground )
                end do
            !if not multiply by the FC factors
            else
                do vib = 0, vibmax
                do vib2 = 0, vibmax
                    hx = nx_1p( lxyz, vib )
                    if ( hx == empty ) cycle
                    hx2 = nx_1p( lxyz2, vib2 )
                    if ( hx2 == empty ) cycle
                    cohfxn( rx,ry, rz ) = cohfxn( rx, ry, rz ) +                 &
                                          h( hx, ground ) * h( hx2, ground ) *   &
                                          fc_gf( 0, vib ) * fc_gf( 0, vib2 )
                end do
                end do
            end if
        end do
        end do
        end do
        end if
        
        !one-particle two-particle
        if ( one_state .and. two_state ) then
        do lx = max( 1-rx, 1 ), min( lattice_x - rx, lattice_x )
        do ly = max( 1-ry, 1 ), min( lattice_y - ry, lattice_y )
        do lz = max( 1-rz, 1 ), min( lattice_z - rz, lattice_z )
            lx2 = lx + rx
            ly2 = ly + ry
            lz2 = lz + rz
            lxyz = nx_lattice( lx, ly, lz )
            lxyz2 = nx_lattice( lx2, ly2, lz2 )
            if ( lxyz == lxyz2 ) cycle            !if on the same molecule, this is zero
            lxyz2v = lxyz                         !free vibration must be on hole site in 1p state
            do vib = 0, vibmax
            do vib2 = 0, vibmax
            do vib2v = 1, vibmax
                hx = nx_1p( lxyz, vib )
                if ( hx == empty ) cycle
                hx2 = nx_2p( lxyz2, vib2, lxyz2v, vib2v )
                if ( hx2 == empty ) cycle
                cohfxn( rx,ry, rz ) = cohfxn( rx, ry, rz ) +                 &
                                      h( hx, ground ) * h( hx2, ground ) *   &
                                      fc_gf( vib2v, vib ) * fc_gf( 0, vib2 ) &
                                      * 2.d0                                !the 2 takes into account 2p-1p interactions (will be the same)
            end do
            end do
            end do
        end do
        end do
        end do
        end if
        
        !two-particle two-particle
        if ( two_state ) then
        do lx = max( 1-rx, 1 ), min( lattice_x - rx, lattice_x )
        do ly = max( 1-ry, 1 ), min( lattice_y - ry, lattice_y )
        do lz = max( 1-rz, 1 ), min( lattice_z - rz, lattice_z )
            lx2 = lx + rx
            ly2 = ly + ry
            lz2 = lz + rz
            lxyz = nx_lattice( lx, ly, lz )
            lxyz2 = nx_lattice( lx2, ly2, lz2 )
        
            !linker type (just like 1p-1p)
            do lxv = 1, lattice_x
            do lyv = 1, lattice_y
            do lzv = 1, lattice_z
                lxyzv = nx_lattice( lxv, lyv, lzv )
                lxyz2v = lxyzv
            do vibv = 1, vibmax
               vib2v = vibv
               if ( lxyz == lxyz2 ) then    !vib must be the same or else fc factor goes to zero
                    do vib = 0, vibmax
                        vib2 = vib
                        hx = nx_2p( lxyz, vib, lxyzv, vibv )
                        if ( hx == empty ) cycle
                        hx2 = nx_2p( lxyz2, vib2, lxyz2v, vib2v )
                        if ( hx2 == empty ) cycle
                        cohfxn( rx,ry, rz ) = cohfxn( rx, ry, rz ) +                 &
                                              h( hx, ground ) * h( hx2, ground )
                    end do
                else                        !vib can be different, now we have fc factors
                    do vib = 0, vibmax
                    do vib2 = 0, vibmax
                        hx = nx_2p( lxyz, vib, lxyzv, vibv )
                        if ( hx == empty ) cycle
                        hx2 = nx_2p( lxyz2, vib2, lxyz2v, vib2v )
                        if ( hx2 == empty ) cycle
                        cohfxn( rx,ry, rz ) = cohfxn( rx, ry, rz ) +                 &
                                              h( hx, ground ) * h( hx2, ground ) *   &
                                              fc_gf( 0, vib ) * fc_gf( 0, vib2 )
                    end do
                    end do
                end if
            end do
            end do
            end do
            end do
            
            !exchange type
            if ( lxyz == lxyz2 ) cycle    !can't be exhange type
            lxyzv = lxyz2
            lxyz2v = lxyz
            do vib = 0, vibmax
            do vib2 = 0, vibmax
            do vibv = 1, vibmax
            do vib2v = 1, vibmax
                    hx = nx_2p( lxyz, vib, lxyzv, vibv )
                    if ( hx == empty ) cycle
                    hx2 = nx_2p( lxyz2, vib2, lxyz2v, vib2v )
                    if ( hx2 == empty ) cycle
                    cohfxn( rx,ry, rz ) = cohfxn( rx, ry, rz ) +                     &
                                          h( hx, ground ) * h( hx2, ground ) *       &
                                          fc_gf( vib2v, vib ) * fc_gf( vibv, vib2 )
            end do
            end do
            end do
            end do
        end do
        end do
        end do
        end if
        
    end do
    end do
    end do
end subroutine
!***************************************************************************************!
!    Write the coherence function
!***************************************************************************************!
subroutine cohfxn_out( cohfxn )
    use exciton_common_local
    implicit none
    
    character(256) f_name
    integer f_no, rx, ry, rz  
    real*8 cohfxn( -lattice_x+1:lattice_x-1, &
                   -lattice_y+1:lattice_y-1, &
                   -lattice_z+1:lattice_z-1 )
    
    !File info
    f_name = trim(task_title)//'_coh.dat'
    f_no = 2
        
    !Open the file and put headers
    open( unit = f_no, file = f_name )
    write(f_no, *) 'Task Title: ', trim(task_title)
    write(f_no, *) 'Polaron Coherence Function'
    write(f_no, *) 'x\y ',(ry, ry = -lattice_y + 1, lattice_y - 1 )
    
    rz = 0
    do rx = -lattice_x+1, lattice_x-1
        write( f_no, * ) rx, ( cohfxn( rx, ry, rz ), ry = -lattice_y + 1, lattice_y - 1 )
    end do
    
    close( f_no )

 End Subroutine