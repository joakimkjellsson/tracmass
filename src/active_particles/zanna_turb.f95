MODULE mod_active_particles

  USE mod_time
  USE mod_traj
  USE mod_grid
  USE mod_param
  USE mod_deformation
  USE mod_vel
  USE mod_loopvars
  
  IMPLICIT none

  REAL                                       :: upr(12,2) = 0
  REAL                                       :: kappa = 1.0
  REAL(DP), ALLOCATABLE, DIMENSION(:,:)      :: mem_param
  REAL(DP), ALLOCATABLE, DIMENSION(:)        :: dt_n, lapu_b, lapu_n, lapv_b, lapv_n
  REAL(DP), DIMENSION(4)                     :: test_sum

CONTAINS

  subroutine active_init
    ALLOCATE ( dt_n(ntracmax) )
    allocate ( lapu_b(ntracmax), lapu_n(ntracmax), lapv_b(ntracmax), lapv_n(ntracmax) )
    ALLOCATE ( mem_param(ntracmax,12) )
    dt_n(:) = 0.
    mem_param(:,:) = 0.
    test_sum(:) = 0.
    lapu_b(:) = 0.
    lapu_n(:) = 0.
    lapv_b(:) = 0.
    lapv_n(:) = 0.
    return
  end subroutine active_init

  subroutine active_ints(ints)
    integer, intent(in) :: ints
    return
  end subroutine active_ints

  subroutine active_ntrac(ntrac)
     
     integer, intent(in) :: ntrac
     integer :: ibp
     real(DP) :: ddx,dlapu,dlapv,lapu1,lapu2,lapv1,lapv2,frac !! temporary variables
     real(DP) :: xx,yy,zz,zuu,zum,zvv,zvm,zx,zy
     integer  :: ic,im,ip,jm,jc,jp,kc
          
     !lapu1 = trajectories(ntrac)%lapu1
     !lapu2 = trajectories(ntrac)%lapu2
     !lapv1 = trajectories(ntrac)%lapv1
     !lapv2 = trajectories(ntrac)%lapv2
     !dlapu = lapu2 - lapu1
     !dlapv = lapv2 - lapv1
     !ic     = trajectories(ntrac)%ib         
     !im = ib-1
     !if(ibm == 0) ibm = IMT
     !jc     = trajectories(ntrac)%jb                   
     !kc     = trajectories(ntrac)%kb                 
     !xx     = trajectories(ntrac)%x1
     !yy     = trajectories(ntrac)%y1
     !zz     = trajectories(ntrac)%z1
     !
     !ip = ic+1
     !if (ic == IMT) ip=1
     !frac = kmt(ic,jc) * kmt(im,jc) * kmt(ip,jc) * kmt(ic,jc-1) * kmt(ic,jc+1)
     !if (frac /= 0) frac=1
     !
     !zuu=(intrpg*uflux(ic,jc,kc,nsp)+intrpr*uflux(ic,jc,kc,nsm))*ff
     !zum=(intrpg*uflux(im,jc,kc,nsp)+intrpr*uflux(im,jc,kc,nsm))*ff
     ! 
     !zvv=(intrpg*vflux(ic,jc  ,kc,nsp)+intrpr*vflux(ic,jc  ,kc,nsm))*ff
     !zvm=(intrpg*vflux(ic,jc-1,kc,nsp)+intrpr*vflux(ic,jc-1,kc,nsm))*ff
     !
     !ddx    = 0.5 * (dxv(ic,jc)+dxv(ic,jc-1)) + 0.5 * (dyu(ic,jc)+dyu(im,jc)) 
     !     
     !! save param for each particle
     !mem_param(ntrac,1) = kappa * ddx**2 * dlapu / dtmin
     !mem_param(ntrac,2) = kappa * ddx**2 * dlapu / dtmin
     !mem_param(ntrac,3) = kappa * ddx**2 * dlapv / dtmin
     !mem_param(ntrac,4) = kappa * ddx**2 * dlapv / dtmin
     
     !mem_param(ntrac,:) = 0.
     
     im = ib-1
     if (im<0) stop 1234
     jm = jb-1
     if (jm<0) stop 2345
     
     zx = x1 - dble(im) ! between 0 and 1                                                                 
     zy = y1 - dble(jm) 
     
     !lapu1 = (intrpg * lapu(im, jb ,kb,nsp) + (1.d0-intrpg) * lapu(im, jb ,kb,nsm)) * (1.d0 - zx) + &
     !        (intrpg * lapu(ib ,jb ,kb,nsp) + (1.d0-intrpg) * lapu(ib ,jb ,kb,nsm)) * zx
     !lapv1 = (intrpg * lapv(ib ,jm ,kb,nsp) + (1.d0-intrpg) * lapv(ib ,jm ,kb,nsm)) * (1.d0 - zy) + &
     !        (intrpg * lapv(ib ,jb ,kb,nsp) + (1.d0-intrpg) * lapv(ib ,jb ,kb,nsm)) * zy
     !print*,intrpg,lapu(im, jb ,kb,nsp),(1.d0-intrpg),lapu(im, jb ,kb,nsm),(1.d0 - zx),&
     !       intrpg,lapu(ib ,jb ,kb,nsp),(1.d0-intrpg),lapu(ib ,jb ,kb,nsm),zx
     !! 
     !! We take laplacian of u at time step nsm (where the particle is currently)
     !! We cant use intrpr and intrpg since they havent been set for this time step yet
     !!
     lapu1 = (0.d0 * lapu(im, jb ,kb,nsp) + (1.d0-0.d0) * lapu(im, jb ,kb,nsm)) * (1.d0 - zx) + &
             (0.d0 * lapu(ib ,jb ,kb,nsp) + (1.d0-0.d0) * lapu(ib ,jb ,kb,nsm)) * zx
     lapv1 = (0.d0 * lapv(ib ,jm ,kb,nsp) + (1.d0-0.d0) * lapv(ib ,jm ,kb,nsm)) * (1.d0 - zy) + &
             (0.d0 * lapv(ib ,jb ,kb,nsp) + (1.d0-0.d0) * lapv(ib ,jb ,kb,nsm)) * zy
     !print*,0.d0,lapu(im, jb ,kb,nsp),(1.d0-0.d0),lapu(im, jb ,kb,nsm),(1.d0 - zx),&
     !       0.d0,lapu(ib ,jb ,kb,nsp),(1.d0-0.d0),lapu(ib ,jb ,kb,nsm),zx
     lapu_b(ntrac) = lapu_n(ntrac)
     lapu_n(ntrac) = lapu1
     lapv_b(ntrac) = lapv_n(ntrac)
     lapv_n(ntrac) = lapv1
     
     !print*,'lapu',lapu1,trajectories(ntrac)%lapu2,trajectories(ntrac)%lapu1
     !print*,'lapu_b, lapu_n ',lapu_b(ntrac),lapu_n(ntrac)
     trajectories(ntrac)%lapu1   = trajectories(ntrac)%lapu2
     trajectories(ntrac)%lapv1   = trajectories(ntrac)%lapv2
     trajectories(ntrac)%lapu2   = lapu1 
     trajectories(ntrac)%lapv2   = lapv2 
          
    return
  end subroutine active_ntrac
  
  subroutine active_niter
     !! Calculate new turbulent fluxes
     !! for each step a particle takes
     !!
     !! RE tensor parameterisation says, after some assumptions, 
     !! du/dt = kappa * dx**2 * d/dt(laplacian(u))
     !!
     !! We then write
     !! u2 = u1 + kappa * dx**2 * dt * d/dt(lapu)
     !! where lapu = d2u/dx2 + d2u/dy2
     !!
     !! Now, the dt * d/dt do not necessarily cancel. 
     !! dt is the time it takes for the particle to 
     !! go travel in this niter iteration. 
     !! d/dt is the inverse of the time between lapu calculations
     !! which maybe between GCM steps. 
     !!
     !! Also, dt is not known yet, since we have not called cross yet.
     !!
     integer  :: ic,im,ip,jm,jc,jp !! temporary indices
     real(DP) :: zup, zvp, zx, zy, zluim, zlui, zdlu, zlvjm, zlvj, zdlv, ddx, zuu, zum, zvv, zvm, zdt, zludu, zlvdv  !! temporary variables  
     real(dp) :: zfrac1, zfrac2, zfrac3, zfrac4
     
     !!
     !! First calculate the turbulent velocity increment
     !!
          
     !! find square of grid box size                                                                                                             
     ddx    = 0.5 * (dxv(ia,ja)+dxv(ia,ja-1)) + 0.5 * (dyu(ia,ja)+dyu(iam,ja))
     
     !! d laplacian(u)
     !zdlu  = trajectories(ntrac)%lapu2 - trajectories(ntrac)%lapu1
     !zdlv  = trajectories(ntrac)%lapv2 - trajectories(ntrac)%lapv1
     !print*,'lap ',zdlu,trajectories(ntrac)%lapu2,trajectories(ntrac)%lapu1
     zdlu = lapu_n(ntrac) - lapu_b(ntrac)
     zdlv = lapv_n(ntrac) - lapv_b(ntrac)
     !print*,'lapu_n, lapu_b, dlapu ',lapu_n(ntrac),lapu_b(ntrac),zdlu
     
#ifdef add_lapu_adv
     !! Add lap(u) * nabla(u)                                                                                                                    
     !! i.e. (d2/dx2 + d2/dy2) u * du/dx + (d2/dx2 + d2/dy2) v * du/dy                                                                           
     !! and  (d2/dx2 + d2/dy2) u * dv/dx + (d2/dx2 + d2/dy2) v * dv/dy                                                                           
     ip = ia+1
     jp = ja+1
     jm = ja-1
     if (ip > imt) ip=ip-imt
     if (jp > jmt) jp=ja
     if (jm < 1) jm=ja
     !zludu = trajectories(ntrac)%lapu1 * (uflux(ip,ja,ka,1) - uflux(iam,ja,ka,1)) / (2.0 * dxv(ia,ja)) + &
     !      & trajectories(ntrac)%lapv1 * (uflux(ia,jp,ka,1) - uflux(ia ,jm,ka,1)) / (2.0 * dyu(ia,ja))
     !zlvdv = trajectories(ntrac)%lapu1 * (vflux(ip,ja,ka,1) - vflux(iam,ja,ka,1)) / (2.0 * dxv(ia,ja)) + &
     !      & trajectories(ntrac)%lapv1 * (vflux(ia,jp,ka,1) - vflux(ia ,jm,ka,1)) / (2.0 * dyu(ia,ja))
     zludu = lapu(ia,ja,ka,1) * (uflux(ip,ja,ka,1) - uflux(iam,ja,ka,1)) / (2.0 * dxv(ia,ja)) + &
           & lapv(ia,ja,ka,1) * (uflux(ia,jp,ka,1) - uflux(ia ,jm,ka,1)) / (2.0 * dyu(ia,ja))
     zlvdv = lapu(ia,ja,ka,1) * (vflux(ip,ja,ka,1) - vflux(iam,ja,ka,1)) / (2.0 * dxv(ia,ja)) + &
           & lapv(ia,ja,ka,1) * (vflux(ia,jp,ka,1) - vflux(ia ,jm,ka,1)) / (2.0 * dyu(ia,ja))
     zludu = zludu * tseas / (dyu(ia,ja) * dzt(ia,ja,ka,1))
     zlvdv = zlvdv * tseas / (dxv(ia,ja) * dzt(ia,ja,ka,1))
#else
     zludu = 0.
     zlvdv = 0.
#endif
     
     !! delta t is length of time step between new velocity fields
     !! however, it should be the time step of the particle, ds 
     !! but we havent called cross and calculated ds yet
     zdt = dt_n(ntrac)/tseas
     
     !if (ia==429 .and. ja==153 .and. ka==64) then
     if (ints == 146027) then
        test_sum(1) = test_sum(1) + zdlu * zdt
        test_sum(2) = test_sum(2) + zdlv * zdt
        !test_sum(3) = test_sum(3) + zdt
        test_sum(4) = test_sum(4) + zdlu
     end if
     !if (y0==152.0) test_sum(4) = -zdlu
     !if (x0==428.0) test_sum(4) = test_sum(4) + zdlu
     test_sum(3) = test_sum(3) + zdt
     !print*,'x0,y0',x0,y0
     !print*,'sum ',test_sum
     
     !! store turbulent velocity increments in upr
     upr(1,:) = kappa * ddx**2 * (zdlu + zludu) * zdt  !! param at ia                                                                            
     upr(2,:) = kappa * ddx**2 * (zdlu + zludu) * zdt  !! param at iam                                                                           
     upr(3,:) = kappa * ddx**2 * (zdlv + zlvdv) * zdt  !! param at ja                                                                            
     upr(4,:) = kappa * ddx**2 * (zdlv + zlvdv) * zdt  !! param at ja-1                                                                          
     upr(5,:) = 0.
     upr(6,:) = 0.
     upr(7,:) = kappa * ddx**2 * (zdlu + zludu) * zdt
     upr(8,:) = kappa * ddx**2 * (zdlu + zludu) * zdt
     upr(9,:) = kappa * ddx**2 * (zdlv + zlvdv) * zdt
     upr(10,:) = kappa * ddx**2 *(zdlv + zlvdv) * zdt
     upr(11,:) = 0.
     upr(12,:) = 0.     
     
     !! integral of all turbulent velocities
     mem_param(ntrac,1) = mem_param(ntrac,1) + upr(1,1)
     mem_param(ntrac,2) = mem_param(ntrac,2) + upr(2,1)
     mem_param(ntrac,3) = mem_param(ntrac,3) + upr(3,1)
     mem_param(ntrac,4) = mem_param(ntrac,4) + upr(4,1)
     
     upr(1,:) = mem_param(ntrac,1)
     upr(2,:) = mem_param(ntrac,2)
     upr(3,:) = mem_param(ntrac,3)
     upr(4,:) = mem_param(ntrac,4)
     upr(7,:) = mem_param(ntrac,1)
     upr(8,:) = mem_param(ntrac,2)
     upr(9,:) = mem_param(ntrac,3)
     upr(10,:) = mem_param(ntrac,4)
     
     !print*,'upr ',upr(1,1),upr(2,1),upr(3,1),upr(4,1)
     
     !! Calculate fluxes on grid cell walls, 
     !! we may need it for later... 
     zuu = (intrpg*uflux(ia ,ja  ,ka,nsp)+intrpr*uflux(ia ,ja  ,ka,nsm))*ff
     zum = (intrpg*uflux(iam,ja  ,ka,nsp)+intrpr*uflux(iam,ja  ,ka,nsm))*ff
     zvv = (intrpg*vflux(ia ,ja  ,ka,nsp)+intrpr*vflux(ia ,ja  ,ka,nsm))*ff
     zvm = (intrpg*vflux(ia ,ja-1,ka,nsp)+intrpr*vflux(ia ,ja-1,ka,nsm))*ff

#ifdef turb_v2
     !! This part checks the magnitude of 
     !! the parameterisation and adjusts 
     !! the ia,ja accordingly. 
     !! This is a pretty disruptive thing to do
     !! Not sure its a good idea. 
     if(zum /= 0.d0) then
        zup = upr(1,1)
     else
        zup = 0.d0
     end if
     if (x0 == dble(iam) .and. zum+zup < 0.d0) then
        ia = iam
        iam = ia-1
        if (iam < 1) iam=iam+imt
        print*,'new ia - ',x0,iam,ia,zum,zup
     end if
     
     if(zuu /= 0.d0) then
        zup = upr(1,1)
     else
        zup = 0.d0
     end if
     if (x0 == dble(ia) .and. zuu+zup > 0.d0) then
        ia = ia+1
        iam = ia-1
        if (ia > imt) ia=ia-imt
        if (iam < 1) iam=iam+imt
        print*,'new ia + ',x0,ia,zuu,zup
     end if
     
     if(zvm /= 0.d0) then 
        zvp = upr(3,1)
     else 
        zvp = 0.d0
     end if
     if (y0 == dble(ja-1) .and. zvm+zvp < 0.d0) then
        ja = ja-1
        print*,'new ja - ',ja,ja-1,zvm,zvp
     end if
     
     if(zvv /= 0.d0) then 
        zvp = upr(3,1)
     else
        zvp = 0.d0
     end if
     if (y0 == dble(ja) .and. zvv+zvp > 0.d0) then
        ja = ja+1
        print*,'new ja + ',ja,zvv,zvp
     end if
#endif
       
!       !! find square of grid box size
!       ddx    = 0.5 * (dxv(ia,ja)+dxv(ia,ja-1)) + 0.5 * (dyu(ia,ja)+dyu(iam,ja))
!       
!       zx     = dmod(x0,1.d0)  !! how far are we from iam
!       zy     = dmod(y0,1.d0)
!       
!       !! store the previous step before calculating a new                              
!       !trajectories(ntrac)%lapu1 = trajectories(ntrac)%lapu2
!       !trajectories(ntrac)%lapv1 = trajectories(ntrac)%lapv2
!       !! find laplacian u at i-1 for current time                                      
!       !zluim = intrpr * lapu(iam,ja  ,ka,nsm) + intrpg * lapu(iam,ja  ,ka,nsp)
!       !zlvjm = intrpr * lapv(ia ,ja-1,ka,nsm) + intrpg * lapv(ia ,ja-1,ka,nsp) 
!       !! find laplacian u at i for current time                                
!       !zlui  = intrpr * lapu(ia,ja,ka,nsm) + intrpg * lapu(ia,ja,ka,nsp)
!       !zlvj  = intrpr * lapv(ia,ja,ka,nsm) + intrpg * lapv(ia,ja,ka,nsp) 
!       !! find laplcian u at zonal position                                             
!       !trajectories(ntrac)%lapu2 = zluim + zx * (zlui - zluim)
!       !trajectories(ntrac)%lapv2 = zlvjm + zy * (zlvj - zlvjm)
!       !! difference from last time                                                       
!       zdlu  = trajectories(ntrac)%lapu2 - trajectories(ntrac)%lapu1
!       zdlv  = trajectories(ntrac)%lapv2 - trajectories(ntrac)%lapv1
!       !! store fluxes
!       !trajectories(ntrac)%dlapu = kappa * ddx**2 * zdlu
!       !trajectories(ntrac)%dlapv = kappa * ddx**2 * zdlv
!       
!#ifdef add_lapu_adv
!       !! Add lap(u) * nabla(u)
!       !! i.e. (d2/dx2 + d2/dy2) u * du/dx + (d2/dx2 + d2/dy2) v * du/dy 
!       !! and  (d2/dx2 + d2/dy2) u * dv/dx + (d2/dx2 + d2/dy2) v * dv/dy
!       ip = ia+1
!       jp = ja+1
!       jm = ja-1
!       if (ip > imt) ip=ip-imt
!       if (jp > jmt) jp=ja
!       if (jm < 1) jm=ja
!       zludu = trajectories(ntrac)%lapu1 * (uflux(ip,ja,ka,1) - uflux(iam,ja,ka,1)) / (2.0 * dxv(ia,ja)) + &
!             & trajectories(ntrac)%lapv1 * (uflux(ia,jp,ka,1) - uflux(ia ,jm,ka,1)) / (2.0 * dyu(ia,ja))
!       zlvdv = trajectories(ntrac)%lapu1 * (vflux(ip,ja,ka,1) - vflux(iam,ja,ka,1)) / (2.0 * dxv(ia,ja)) + &
!             & trajectories(ntrac)%lapv1 * (vflux(ia,jp,ka,1) - vflux(ia ,jm,ka,1)) / (2.0 * dyu(ia,ja))
!#else
!       zludu = 0.
!       zlvdu = 0.
!#endif              
!       
!       !lapu1  = trajectories(ntrac)%lapu1                                                                                   
!       !lapu2  = trajectories(ntrac)%lapu2                                                                             
!       !! dlapu is the difference in laplacian of u                                                                          
!       !! from the beginning of last gcm time step                                                                           
!       !! to the current gcm time step                                                                                      
!       !! If a particle stays within the same grid cell                                                                      
!       !! dlapu = 0                                                                                                           
!       !zdlu  = trajectories(ntrac)%dlapu / (ngcm * 3600.)                                         
!       !zdlv  = trajectories(ntrac)%dlapv / (ngcm * 3600.)                                                          
!       
!       !! add same fluxes to u and v
!       !! so that we dont change the divergence
!       !! We could calculate div and add as a param for w 
!       zdt = dtmin/tseas
!       
!       upr(1,:) = kappa * ddx**2 * (zdlu + zludu) * zdt  !! param at ia
!       upr(2,:) = kappa * ddx**2 * (zdlu + zludu) * zdt  !! param at iam
!       upr(3,:) = kappa * ddx**2 * (zdlv + zlvdv) * zdt  !! param at ja
!       upr(4,:) = kappa * ddx**2 * (zdlv + zlvdv) * zdt  !! param at ja-1
!       upr(5,:) = 0.
!       upr(6,:) = 0.
!       upr(7,:) = kappa * ddx**2 * (zdlu + zludu) * zdt 
!       upr(8,:) = kappa * ddx**2 * (zdlu + zludu) * zdt 
!       upr(9,:) = kappa * ddx**2 * (zdlv + zlvdv) * zdt 
!       upr(10,:) = kappa * ddx**2 *(zdlv + zlvdv) * zdt 
!       upr(11,:) = 0.
!       upr(12,:) = 0.
       
       !! Adjust upr so that we dont change the sign of 
       !! uflux and vflux on grid box walls
       !zfrac1 = min( abs(upr(1,1)), abs(0.99*zuu/upr(1,1)) )
       !zfrac2 = min( abs(upr(2,1)), abs(0.99*zum/upr(2,1)) )
       !zfrac3 = min( abs(upr(3,1)), abs(0.99*zvv/upr(3,1)) )
       !zfrac4 = min( abs(upr(4,1)), abs(0.99*zvm/upr(4,1)) )
       
       if (upr(1,1) /= 0.d0) then
          if (abs(upr(1,1)) > abs(zuu)) then
             zfrac1 = 0.99 * abs(zuu/upr(1,1))
          else
             zfrac1 = 1.d0
          end if
       else
          zfrac1 = 0.d0
       end if
       
       if (upr(2,1) /= 0.d0) then
          if (abs(upr(2,1)) > abs(zum)) then
             zfrac2 = 0.99 * abs(zum/upr(2,1))
          else
             zfrac2 = 1.d0
          end if
       else
          zfrac2 = 0.d0
       end if
       
       if (upr(3,1) /= 0.d0) then
          if (abs(upr(3,1)) > abs(zvv)) then
             zfrac3 = 0.99 * abs(zvv/upr(3,1))
          else
             zfrac3 = 1.d0
          end if
       else
          zfrac3 = 0.d0
       end if
       
       if (upr(4,1) /= 0.d0) then
          if (abs(upr(4,1)) > abs(zvm)) then
             zfrac4 = 0.99 * abs(zvm/upr(4,1))
          else
             zfrac4 = 1.d0
          end if
       else
          zfrac4 = 0.d0
       end if 
       
       if (ntrac==27628) print*,' uu ',zuu,zum,zvv,zvm
       if (ntrac==27628) print*, 'upr ',upr(1,1),upr(2,1),upr(3,1),upr(4,1)
       
       !upr(1,:) = upr(1,:) * zfrac1 !! param at ia                                                                                  
       !upr(2,:) = upr(2,:) * zfrac2 !! param at iam                                                                                 
       !upr(3,:) = upr(3,:) * zfrac3 !! param at ja                                                                                  
       !upr(4,:) = upr(4,:) * zfrac4 !! param at ja-1                                                                                
       !! Make sure to add same fluxes on i and i-1, and j and j-1
       !! so that divergence in grid cells remain unchanged
       if (abs(zfrac1*upr(1,1)) <= abs(zfrac2*upr(2,1))) then
          upr(1,:) = zfrac1*upr(1,:) 
          upr(2,:) = upr(1,:)
       else
          upr(1,:) = zfrac2*upr(2,:)
          upr(2,:) = zfrac2*upr(2,:)
       end if
       if (abs(zfrac3*upr(3,1)) <= abs(zfrac4*upr(4,1))) then
          upr(3,:) = zfrac3*upr(3,:) 
          upr(4,:) = upr(3,:)
       else
          upr(3,:) = zfrac4*upr(4,:)
          upr(4,:) = zfrac4*upr(4,:)
       end if
       upr(5,:) = 0. 
       upr(6,:) = 0.
       !upr(7,:) = upr(7,:) * zfrac1
       !upr(8,:) = upr(8,:) * zfrac2
       !upr(9,:) = upr(9,:) * zfrac3
       !upr(10,:) = upr(10,:) * zfrac4
       upr(7,:) = upr(1,:) 
       upr(8,:) = upr(2,:)
       upr(9,:) = upr(3,:) 
       upr(10,:) = upr(4,:)
       upr(11,:) = 0.
       upr(12,:) = 0.
       
       !! set turbulent velocities to zero
       !! if original velocity is zero
       if (zuu == 0.d0 .or. zum == 0.d0) then
          upr(1,:) = 0.
          upr(2,:) = 0.
          upr(7,:) = 0.
          upr(8,:) = 0.
       end if
       if (zvv == 0.d0 .or. zvm == 0.d0) then
          upr(3,:) = 0.
          upr(4,:) = 0.
          upr(9,:) = 0.
          upr(10,:) = 0.
       end if
       if (ntrac==27628) print*,'upr 2 ',upr(1,1),upr(2,1),upr(3,1),upr(4,1)

    return
  end subroutine active_niter
  
   subroutine active_niter_2
      dt_n(ntrac) = ds * dxyz
   return
   end subroutine active_niter_2

END MODULE mod_active_particles



