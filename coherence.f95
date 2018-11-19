!***************************************************************************************!
!                                    This subroutine calculates the coherence function
!***************************************************************************************!
subroutine cms_cohfxn( cohfxn )
    use common_variables
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
     if (eval(ground)> -3 .and. eval(ground)< -2) then
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
            lxyz2v = lxyz                               !free vibration must be on hole site in 1p state
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
end if
end subroutine
!***************************************************************************************!
!    Write the coherence function
!***************************************************************************************!
subroutine cohfxn_out( cohfxn )
     use common_variables
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