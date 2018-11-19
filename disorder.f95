
!***************************************************
!***************************************************
!***************************************************
                !  This subroutine adds diagonal disorder 
!***************************************************
!***************************************************
!***************************************************
 subroutine short_short(config)
 use common_variables
 implicit none

integer vx,vx1,vy,vy1,vz,vz1,v1,v2,vxyz,vxyz1,hx1,hx2,config


do vx = 1, lattice_x    
        do vy = 1, lattice_y
        do vz = 1, lattice_z
           do v1=0,vibmax  
          vxyz = nx_lattice( vx, vy, vz )
            hx1= nx_1p(vxyz,v1) 
            if (hx1==empty) cycle
             h(hx1,hx1) = h(hx1,hx1) + disorder (config,vx,vy)
             write(7,*) 'h(',hx1,hx1,')', h(hx1,hx1)

 do vx1=1,lattice_x
            do vy1=1,lattice_y
            do vz1=1,lattice_z
            do v2=1,vibmax
             vxyz1=nx_lattice(vx1,vy1,vz1)
          hx2=nx_2p(vxyz,v1,vxyz1,v2)
          if (hx2==empty) cycle
          
          
          h(hx2,hx2) = h(hx2,hx2) + disorder (config,vx,vy)
          write(7,*) 'h(',hx2,hx2,')', h(hx2,hx2)
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        end do     
      !  end if

       
       
 end subroutine


!***************************************************************************************!
!    This subroutine introduces disorder 
!***************************************************************************************!

        
        subroutine disorder_elements()
        use common_variables
        implicit none
        
        integer vx,vy
         real *8 rand
real*8, external :: dlarnd
integer config, idist, i, j
    parameter ( idist = 3 )                      
    integer :: iseed(4) = (/47, 3093,1041,77/)   
        
         if ( .not. allocated( disorder ) ) then 
            allocate( disorder ( conf_max, lattice_x, lattice_y ))
    end if
       
       do config =1, conf_max
       
        do vx=1,lattice_x
      !  rand = dlarnd( idist, iseed ) 
      !  rand = rand * 2.01625
        do vy=1,lattice_y
        rand = dlarnd( idist, iseed ) 
        rand = rand * 1.152142857
        !rand = rand * 1.152247143
        disorder (config,vx,vy) = rand
        write (10,*) 'disorder (', config , vx , vy , ')', rand
        end do
end do
end do

end subroutine    

!*************************************************************************************************!
!    This subroutine introduces a dopant anion which interacts coulombically with the hole in the polymer aggregate
!*************************************************************************************************!



subroutine doping_disorder()
 use common_variables
        implicit none
        
        integer vx,vy,dx,dy,config

        if ( .not. allocated( D ) ) then 
            allocate( D (conf_max, lattice_x, lattice_y ))
    end if
    
        if ( .not. allocated( Force ) ) then 
            allocate( Force ( conf_max, lattice_x, lattice_y ))
    end if
!    Force = 0
       do config =1, conf_max
       
      !  D(config,3,3) = r 
        do vx=1,lattice_x
        do vy=1,lattice_y
      
       dx = abs(vx-3)
       dy = abs(vy-3)
       
        D(config, vx,vy) =  sqrt((dintra*dy)**2 + (dinter*dx)**2 + r**2)
       Force(config,vx,vy) = (-8.3)/(D(config,vx,vy))
        print*, 'D (', config , vx , vy , ')', D(config,vx,vy)
    !    print*, 'Force (', config , vx , vy , ')', Force(config,vx,vy)
        end do
        end do
        end do
        end subroutine
       
  subroutine paracrystalline_SS()
        use common_variables
        implicit none
        
         integer lx1,lx2, ly1,ly2, dx, dy
       
         real *8 rand,g
        real *8 :: sum_te = 0
        
real*8, external :: dlarnd
integer config, idist, i, j
    parameter ( idist = 3 )  
                
    integer :: iseed(4) = (/27, 3091,1043,67/)   
  if ( .not. allocated( pdisorder ) ) then 
            allocate( pdisorder ( conf_max,lattice_x,lattice_y, lattice_x, lattice_y))
    end if
    rand = 0
  
     
  do config=1, conf_max

 !rand = dlarnd( idist, iseed )
! rand = rand* 0.247 
! g = 0.247/3.8    
  do lx1=1, lattice_x
  do ly1=1, lattice_y-1
  
  

 ly2=ly1+1
  rand = dlarnd( idist, iseed )
  rand = rand* 0.076
  g = 0.076/3.8    

 
 
 pdisorder (config,lx1,ly1,lx1, ly2) = rand
 pdisorder (config,lx1,ly2,lx1, ly1) = rand
 
 
  
!  print*, 'pdisorder (', config,lx1,ly1,lx1, ly2, ')', pdisorder (config,lx1,ly1,lx1, ly2)
 ! print*, 'pdisorder (', config,lx1,ly2,lx1, ly1, ')', pdisorder (config,lx1,ly2,lx1, ly1)

   
     ! ly1= ly1+ 1
     ! ly2 = ly2 +1  
    end do
   end do
     end do 
      
        end subroutine