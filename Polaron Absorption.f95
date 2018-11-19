!***************************************************************************************!
!    Calculate the CMS oscillator strength.
!***************************************************************************************!
subroutine cms_osc()
    use common_variables
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
          if (eval(ground)> -3 .and. eval(ground)< -2) then
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
      end if  
    end do
end subroutine


!***************************************************************************************!
!    Calculate the CMS oscillator strength sums
!***************************************************************************************!
subroutine cms_osc_sum( sumx, sumy )
    use common_variables
    implicit none
    
    real*8, intent(inout) :: sumx, sumy
    integer excited, ground
    parameter ( ground = 1 )
    real*8 trans_e
    
    do excited = 2, kount2
      if (eval(ground)> -3 .and. eval(ground)< -2) then
        trans_e = eval( excited ) - eval( ground )
        sumx = sumx + trans_e * osc_x( excited )
        sumy = sumy + trans_e * osc_y( excited )   
        end if
    end do

end subroutine
!***************************************************************************************!
!    Calculate the CMS spectrum
!***************************************************************************************!
subroutine cms_spec( ab_x, ab_y )
    use common_variables
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
          if (eval(ground)> -3 .and. eval(ground)< -2) then
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
 end if
        end do        
    end do
  
end subroutine
!***************************************************************************************!
!    Write the CMS spectrum
!***************************************************************************************!
subroutine cms_out( ab_x, ab_y )
    use common_variables
    implicit none
    
    character(256) f_name
    integer f_no, spec_point
    real*8 spec_start, spec_end, energy,                                 &
           ab_x( spec_step ), ab_y( spec_step )
        
    !File info
    f_name = trim(task_title)//'_cms1.dat'
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
        write( f_no, '(5f14.7)' ) energy,                                           &
                                  (ab_x( spec_point ) + ab_y( spec_point )), &
                                  ab_x( spec_point ) , ab_y( spec_point )
    end do
    close( f_no )

 End Subroutine