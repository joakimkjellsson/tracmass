SUBROUTINE readfields
  
   USE mod_precdef
   USE mod_param
   USE mod_vel
   
   USE mod_time
   USE mod_grid
   USE mod_name
   USE mod_vel
   USE mod_traj
   USE mod_getfile
   use mod_seed
   use mod_tempsalt
   USE mod_deformation
   USE mod_laplacian
   
#ifdef tempsalt
   USE mod_dens
   USE mod_stat
#endif
   IMPLICIT none
   ! ==========================================================================
   ! === Read velocity, temperature and salinity for INALT60.L120 configuration ===
   ! ==========================================================================
   ! Run each time step
   ! --------------------------------------------------------------------------
   ! The following arrays will be populated:
   !
   ! uflux    - Zonal volume flux (U point)
   ! vflux    - Meridional volume flux (V point)
   !
   ! If run with tempsalt option, the following are also set
   ! tem      - Temperature (T point) 
   ! sal      - Salinity (T point)
   ! rho      - Potential density (T point)
   !
   ! --------------------------------------------------------------------------
   
   ! = Loop variables
   INTEGER                                       :: i, j, k ,kk, im, ip, jm, jp, imm, ii, jmm, jpp, l
   INTEGER                                       :: kbot,ktop 
   INTEGER, SAVE                                 :: itime, fieldStep
   INTEGER, SAVE                                 :: ntempus=0,ntempusb=0,nread
   ! = Variables used for getfield procedures
   CHARACTER (len=200)                           :: fieldFile, medfieldFile, timestamp, tmpstr
   ! = Variables for filename generation
   CHARACTER (len=200)                           :: dataprefix
   REAL(DP), ALLOCATABLE, DIMENSION(:,:)         :: zstot,zstou,zstov,abyst,abysu,abysv
   REAL(DP), ALLOCATABLE, DIMENSION(:,:,:)       :: xxx
   REAL*4                                      :: dd,hu,hv,uint,vint,zint,hh,h0
   REAL*4, ALLOCATABLE, DIMENSION(:,:,:),SAVE    :: u_m, v_m
   
#ifdef tempsalt
   REAL*4, ALLOCATABLE, DIMENSION(:)           :: rhozvec, depthzvec, latvec
   REAL*4, ALLOCATABLE, DIMENSION(:)           :: tmpzvec, salzvec
#endif
   
   INTEGER, ALLOCATABLE, DIMENSION(:,:),SAVE          :: fileMon, fileDay
   CHARACTER (len=200), ALLOCATABLE, DIMENSION(:),SAVE :: file_timestamp
    
   LOGICAL                                       :: around
 
!---------------------------------------------------------------
   
   !
   ! Allocate variables 
   !
   alloCondUVW: if(.not. allocated (zstot)) then
      allocate ( zstot(imt,jmt),zstou(imt,jmt),zstov(imt,jmt) )
      allocate ( xxx(imt,jmt,km))
#ifdef tempsalt
      allocate ( tmpzvec(km), salzvec(km), rhozvec(km), depthzvec(km), latvec(km))
#endif
   endif alloCondUVW
 
   call datasetswap ! Swap between current and previous step
   vort(:,:,:,1) = vort(:,:,:,2)
   hdiv(:,:,:,1) = hdiv(:,:,:,2)
   lapu(:,:,:,1) = lapu(:,:,:,2)
   lapv(:,:,:,1) = lapv(:,:,:,2)
   call updateClock 
 
! === Initialising fields ===
   if (ints == intstart) then      
      
      ! INALT60 data with 4h frequency
      allocate( file_timestamp(16), fileDay(16,2), fileMon(16,2) )
      fileDay(1:16,1) = (/  1,  6, 31, 25, 22, 16, 11,  5, 30, 25, 19, 13,  8,  2, 27, 22 /)
      fileDay(1:16,2) = (/  5, 30, 24, 21, 15, 10,  4, 29, 24, 18, 12,  7,  1, 26, 21, 31 /)
      fileMon(1:16,1) = (/  1,  1,  1,  2,  3,  4,  5,  6,  6,  7,  8,  9, 10, 11, 11, 12 /)
      fileMon(1:16,2) = (/  1,  1,  2,  3,  4,  5,  6,  6,  7,  8,  9, 10, 11, 11, 12, 12 /)
      
      tmpstr = "xxxxxxxx_xxxxxxxx"
      !print*,'tmpstr: ',tmpstr
      do ii=1,16
         write(tmpstr(5:8 ),'(i2.2,i2.2)') fileMon(ii,1),fileDay(ii,1)
         write(tmpstr(14:17),'(i2.2,i2.2)') fileMon(ii,2),fileDay(ii,2)
         !print*,'tmpstr: ',trim(tmpstr)
         file_timestamp(ii) = trim(tmpstr)
         !print*,'file_timestamp(ii): ',trim(file_timestamp(ii))
      end do      
   
      !  
      ! Find the files for the current step  
      ! 
      fieldStep = 0
      itime = 1
      !print*,'currTime,fileTime,diff: ',currMon*100 + currDay,fileMon(itime,2)*100 + fileDay(itime,2),&
      !& currMon*100 + currDay - (fileMon(itime,2)*100 + fileDay(itime,2))
      do while (currMon*100 + currDay - (fileMon(itime,2)*100 + fileDay(itime,2)) > 0)
         !print*,'currTime,fileTime,diff,itime: ',currMon*100 + currDay,fileMon(itime,2)*100 + fileDay(itime,2),&
         !& currMon*100 + currDay - (fileMon(itime,2)*100 + fileDay(itime,2)),itime
         itime = itime + 1
      end do            
             
   end if
   
   fieldStep = fieldStep + 1   
   print*,'fieldStep: ',fieldStep
   print*,'fileMon(itime,1)*100+fileDay(itime,1)',fileMon(itime,1)*100+fileDay(itime,1)
   print*,'hej'
   if ( (fieldStep > 150 .and. ( fileMon(itime,1)*100+fileDay(itime,1) > 101 .and. &
                               & fileMon(itime,1)*100+fileDay(itime,1) < 1222 ) ) .or. & 
      & (fieldStep > 30  .and. ( fileMon(itime,1)*100+fileDay(itime,1) == 101)  ) .or. &
      & (fieldStep > 60  .and. ( fileMon(itime,1)*100+fileDay(itime,1) == 1222) ) ) then
      fieldStep = 1
      itime = itime + 1
      if (itime > 16) then
         itime = 1
      end if
   end if 
   
   ncTpos = fieldStep
   
   !print*,'itime:', itime
   timestamp = trim(file_timestamp(itime))
   !print*,'timestamp',timestamp
   write(timestamp(1:4 ),'(i4)') currYear
   write(timestamp(10:13),'(i4)') currYear
   !print*,'timestamp:',timestamp
   dataprefix='2_INALT60.L120-KRS0020_4h_'//trim(timestamp)
   !print*,'dataprefix:',dataprefix
   fieldFile = trim(inDataDir)//trim(dataprefix)
   
   print*,' Date: ',currYear, currMon, currDay
   print*,' Time: ',currHour
   print*,' File: ',fieldFile
   print*,' Step: ',fieldStep
   
   ! Read SSH
   hs(:,     :, nsp) = get2DfieldNC(trim(fieldFile)//'_grid_T.nc', 'sossheig')
   hs(imt+1, :, nsp) = hs(1,:,nsp)
   
   ! Depth at U, V, T points as 2D arrays
   allocate ( abyst(imt, jmt) , abysu(imt, jmt) , abysv(imt, jmt) )
   
   abyst = sum(dzt0(:,:,:), dim=3)
   abysu = sum(dzu(:,:,:,1), dim=3)
   abysv = sum(dzv(:,:,:,1), dim=3)
   
   ! Calculate SSH/depth
   where (abyst /= 0)
      zstot = hs(:imt,:jmt,nsp)/abyst + 1
   elsewhere
      zstot = 0.d0
   end where
   
   where (abysu /= 0)
      zstou = 0.5*(hs(:imt,:jmt,nsp)+hs(2:imt+1,:jmt,nsp))/abysu + 1
   elsewhere
      zstou = 0.d0
   end where
   
   where (abysv /= 0)
      zstov = 0.5*(hs(:imt,:jmt,nsp)+hs(:imt,2:jmt+1,nsp))/abysv + 1
   elsewhere
      zstov = 0.d0
   end where
 
   ! Read temperature 
#if defined tempsalt 
   xxx(:,:,:) = get3DfieldNC(trim(fieldFile)//'_grid_T.nc', 'votemper')
   tem(:,:,:,nsp) = xxx(:,:,km:1:-1)
   
   ! Read salinity
   xxx(:,:,:) = get3DfieldNC(trim(fieldFile)//'_grid_T.nc', 'vosaline')
   sal(:,:,:,nsp) = xxx(:,:,km:1:-1)
   
   ! Calculate potential density
   depthzvec = 0.
   do j=1,jmt
      latvec=-80+1./12.*float(j+subGridJmin-1)
      do i=1,IMT
         tmpzvec = tem(i,j,:,nsp)
         salzvec = sal(i,j,:,nsp)
         call statvd(tmpzvec, salzvec, rhozvec ,km ,depthzvec ,latvec)
         rho(i,j,:,nsp)=rhozvec - 1000.
      end do
   end do
#endif     
   
   ! Read u, v
   uvel = get3DfieldNC(trim(fieldFile)//'_grid_U.nc', 'vozocrtx')
   vvel = get3DfieldNC(trim(fieldFile)//'_grid_V.nc', 'vomecrty')
   
   !
   ! Calculate zonal and meridional volume flux
   !
   ! Weight by (1 + ssh / depth)
   ! This is only an approximation of what NEMO really does
   ! but is accurate within 1% 
   !
   
   do k = 1, km
   do j = 1, jmt
   do i = 1, imt
      dzt(i,j,k,nsp) = dzt0(i,j,k) * zstot(i,j)
   end do
   end do
   end do

#ifdef nomean   
   !! flip u,v upside down
   !! use uflux, vflux as temporary arrays
   where (uvel == 0)
      u_m = 0
   end where
   where (vvel == 0)
      v_m = 0
   end where
   
   uflux(1:imt,1:jmt,1:km,nsp) = uvel(1:imt,1:jmt,1:km) - u_m(1:imt,1:jmt,1:km)
   vflux(1:imt,1:jmt,1:km,nsp) = vvel(1:imt,1:jmt,1:km) - v_m(1:imt,1:jmt,1:km)

#else
   uflux(1:imt,1:jmt,1:km,nsp) = uvel(1:imt,1:jmt,1:km) 
   vflux(1:imt,1:jmt,1:km,nsp) = vvel(1:imt,1:jmt,1:km) 
#endif
      
   uvel(:,:,:) = 0.
   vvel(:,:,:) = 0.
   do k = 1, km
   do j = 1, jmt
   do i = 1, imt
      uvel(i,j,km+1-k) = uflux(i,j,k,nsp) 
      vvel(i,j,km+1-k) = vflux(i,j,k,nsp) 
   enddo
   enddo
   enddo
   
   !! calculate volume fluxes
   uflux(:,:,:,nsp) = 0.
   vflux(:,:,:,nsp) = 0.
   do k = 1, km
   do j = 1, jmt
   do i = 1, imt
      uflux(i,j,km+1-k,nsp) = uvel(i,j,km+1-k) * dyu(i,j) * dzu(i,j,km+1-k,1) * zstou(i,j)
      vflux(i,j,km+1-k,nsp) = vvel(i,j,km+1-k) * dxv(i,j) * dzv(i,j,km+1-k,1) * zstov(i,j)
   enddo
   enddo
   enddo
   
   !! calculate laplacian of u,v
   call laplacian
      
   ! Check that volume fluxes are zero below sea floor
   do i=1,IMT
   do j=1,JMT
   do k=1,KM
   if(k > kmv(i,j) .and. vflux(i,j,km+1-k,nsp) /= 0.) then
      print *,'vflux=',vflux(i,j,km+1-k,nsp),vvel(i,j,k),i,j,k,kmv(i,j),nsp
      stop 4966
   endif
   if(k > kmu(i,j) .and. uflux(i,j,km+1-k,nsp) /= 0.) then
      print *,'uflux=',uflux(i,j,km+1-k,nsp),uvel(i,j,k),i,j,k,kmu(i,j),nsp
      stop 4967
   endif
   enddo
   enddo
   enddo
   
   
#ifdef drifter
   ! average velocity/transport to simulate drifter trajectories
   kbot=65 ; ktop=66 ! number of surface layers to integrate over 
   uint=0. ; vint=0. ; zint=0.
   do k=kbot,ktop
      uint = uint + uflux(:,:,k,nsp) ! integrated transport
      vint = vint + vflux(:,:,k,nsp)
      zint = zint + dz(k)          ! total depth of drougued drifter
   end do
   ! weighted transport for each layer
   do k=kbot,KM
      uflux(:,:,k,nsp) = uint*dz(k)/zint 
      vflux(:,:,k,nsp) = vint*dz(k)/zint
   enddo
#endif

#ifdef initxyt
   ! Set the initial trajectory positions
   !ijkst(:,5)=ntimask(ntempus,:)
#ifdef orca025l75h6
   if( mod(ints,24/ngcm*5) == 1 .or. ints <= 2) ntempus=ntempus+1
   if(ntempus /= ntempusb .and. ntempus <= NTID) then
      ntempusb=ntempus
      !print *,'ints=',ints,' ntempus=',ntempus,' ntempusb=',ntempusb
#else
   if(ints.le.NTID) then
#endif
   
      do ntrac=1,ijkmax
      if(trajinit(ntempus,ntrac,3) /= 0.) then
         ijkst(ntrac,4)=0
         ijkst(ntrac,5)=5
         ijkst(ntrac,6)=ijkst(ntrac,6)+1
         do l=1,3
            ijkst(ntrac,l)=trajinit(ntempus,ntrac,l)+1
            trj(ntrac,l)=trajinit(ntempus,ntrac,l)
            !if(l.eq.1) print *,ntrac,float(ijkst(ntrac,l)-1),trj(ntrac,l),float(ijkst(ntrac,l))
            if(trj(ntrac,l).gt.float(ijkst(ntrac,l)) .or. trj(ntrac,l).lt.float(ijkst(ntrac,l)-1)) then
               print *,l,ntrac,float(ijkst(ntrac,l)-1),trj(ntrac,l),float(ijkst(ntrac,l))
               stop 3946
            endif
         enddo
      else
         ijkst(ntrac,5)=0
         ijkst(ntrac,6)=0
      endif
      enddo
   endif

#ifdef orca025l75h6
#endif
   if( mod(ints,24/ngcm*5).ne.1 .and. ints.gt.2) then
      ijkst=0 
   endif
#endif

   return
   
end subroutine readfields


