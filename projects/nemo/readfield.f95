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
   
   USE mod_dens
   USE mod_stat
   
   IMPLICIT none
   
   ! ==========================================================================
   ! === Read velocity, temperature and salinity for ORCA0083 configuration ===
   ! ==========================================================================
   ! Subroutine to read the ocean state from ORCA0083 config
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
   INTEGER                                       :: ichar, itime, fieldStep
   INTEGER, SAVE                                 :: ntempus=0,ntempusb=0,nread
   ! = Variables used for getfield procedures
   CHARACTER (len=200)                           :: fieldFile, medfieldFile, umFile, vmFile
   CHARACTER (len=200)                           :: tFile, uFile, vFile
   CHARACTER (len=100)                           :: tmpstr
   ! = Variables for filename generation
   CHARACTER (len=200)                           :: dataprefix, timestamp
   REAL(DP), ALLOCATABLE, DIMENSION(:,:)         :: zstot,zstou,zstov,abyst,abysu,abysv
   REAL(DP), ALLOCATABLE, DIMENSION(:,:,:)       :: xxx
   REAL*4                                      :: dd,hu,hv,uint,vint,zint,hh,h0
   REAL*4, ALLOCATABLE, DIMENSION(:,:,:),SAVE    :: u_m, v_m
   
   INTEGER, ALLOCATABLE, DIMENSION(:,:),SAVE          :: fileMon, fileDay
   CHARACTER (len=200), ALLOCATABLE, DIMENSION(:),SAVE :: file_timestamp
   
   INTEGER, PARAMETER                            :: NTID=73
   INTEGER, PARAMETER                            :: IJKMAX2=7392 ! for distmax=0.25 and 32 days

   INTEGER,  SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: ntimask
   REAL*4, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: trajinit
   
   REAL*4, ALLOCATABLE, DIMENSION(:)           :: rhozvec, depthzvec, latvec
   REAL*4, ALLOCATABLE, DIMENSION(:)           :: tmpzvec, salzvec

   LOGICAL                                       :: around, readMean = .false., readTS = .true., &
                                                    vvl = .true., readBio = .false., readSSH = .true.
 
!---------------------------------------------------------------

if (ints == intstart) then
   
   if (readMean) then
      !! Read mean  
      if(.not. allocated (u_m)) then
         allocate ( u_m(imt,jmt,km), v_m(imt,jmt,km) )
         print*,' Read mean fields '
         
         umFile = '/group_workspaces/jasmin2/aopp/joakim/ORCA0083-N001/mean/ORCA0083-N01_1978-2010m00U.nc'
         vmFile = '/group_workspaces/jasmin2/aopp/joakim/ORCA0083-N001/mean/ORCA0083-N01_1978-2010m00V.nc'
         u_m(1:imt,1:jmt,1:km) = get3DfieldNC(umFile,'vozocrtx')
         v_m(1:imt,1:jmt,1:km) = get3DfieldNC(vmFile,'vomecrty')
      end if
   end if
   
      ! INALT60 data with 4h frequency
      allocate( file_timestamp(16), fileDay(16,2), fileMon(16,2) )
      fileDay(1:16,1) = (/  1,  6, 31, 25, 22, 16, 11,  5, 30, 25, 19, 13,  8,  2, 27, 22 /)
      fileDay(1:16,2) = (/  5, 30, 24, 21, 15, 10,  4, 29, 24, 18, 12,  7,  1, 26, 21, 31 /)
      fileMon(1:16,1) = (/  1,  1,  1,  2,  3,  4,  5,  6,  6,  7,  8,  9, 10, 11, 11, 12 /)
      fileMon(1:16,2) = (/  1,  1,  2,  3,  4,  5,  6,  6,  7,  8,  9, 10, 11, 11, 12, 12 /)
      
      tmpstr = "YYYYMMDD_YYYYMMDD"
      print*,'tmpstr: ',tmpstr
      do ii=1,16
         write(tmpstr(5:8)  ,'(i2.2,i2.2)') fileMon(ii,1),fileDay(ii,1)
         write(tmpstr(14:17),'(i2.2,i2.2)') fileMon(ii,2),fileDay(ii,2)
         print*,'tmpstr: ',trim(tmpstr)
         file_timestamp(ii) = trim(tmpstr)
         print*,'file_timestamp(ii): ',trim(file_timestamp(ii))
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

#ifdef initxyt
   ! 
   ! Allocate variables necessary for drifter simulation
   !
   alloCondGrid: if ( .not. allocated (ntimask) ) then
      allocate ( ntimask(NTID,IJKMAX2,3) , trajinit(NTID,IJKMAX2,3) )
   endif alloCondGrid
#endif
   
   !
   ! Allocate variables 
   !
   alloCondUVW: if(.not. allocated (zstot)) then
      allocate ( zstot(imt,jmt),zstou(imt,jmt),zstov(imt,jmt) )
      allocate ( xxx(imt,jmt,km))
      if (readTS) then
         allocate ( tmpzvec(km), salzvec(km), rhozvec(km), depthzvec(km), latvec(km))
      end if
   endif alloCondUVW
 
   call datasetswap ! Swap between current and previous step
   vort(:,:,:,1) = vort(:,:,:,2)
   hdiv(:,:,:,1) = hdiv(:,:,:,2)
   lapu(:,:,:,1) = lapu(:,:,:,2)
   lapv(:,:,:,1) = lapv(:,:,:,2)
   call updateClock 
 
! === Initialising fields ===
   initFieldcond: if(ints == intstart) then
   
#ifdef initxyt
      ! Time for individual start positions
      if(IJKMAX2 == 7392) open(84,file=trim(inDataDir)//'topo/masktime_32_025', &
                               form='unformatted')
      read(84) trajinit
      close(84)
      j=0
      do k=1,NTID
         do i=1,IJKMAX2
            if(trajinit(k,i,3) /= 0.) then
               j=j+1
#if orca025l75h6
               trajinit(k,i,3)=float(kst2)-0.5
               !   print *,j,trajinit(k,i,:)
#endif
            endif
         enddo
      enddo
      ijkst=0
      if(j /= IJKMAX2) then
         stop 4396
      endif
#endif
   
   endif initFieldcond
   
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
   
   !
   ! Find the files for the current step
   !
   
   !physPrefixForm = 'RUNID_TSTSTSTS_GRIDX'
   !RunID = '2_INALT60.L120-KRS0020_4h'
   !tGridName = 'grid_T'
   !fileSuffix = '.nc'
   
   !ssh_name = 'sossheig'
   !ueul_name = 'vozocrtx'
   !veul_name = 'vomecrty'
   
   !physTracerNames(1) = 'votemper'
   !physTracerNames(2) = 'vosaline'
   
   
   print*,physPrefixForm
   print*,file_timestamp,itime
   timestamp = trim(file_timestamp(itime))
   print*,'timestamp',timestamp
   
   ichar = INDEX(physPrefixForm,'RUNID')
   do while (ichar /= 0)
      physPrefixForm = trim(physPrefixForm(:ichar-1))//trim(RunID)//trim(physPrefixForm(ichar+5:))
      print*,physPrefixForm
      ichar = INDEX(physPrefixForm,'RUNID')
   end do
   
   ichar = INDEX(physPrefixForm,'YYYY')
   print*,physPrefixForm,ichar
   do while (ichar /= 0)
      write(physPrefixForm(ichar:ichar+3),'(i4)') currYear
      ichar = INDEX(physPrefixForm,'YYYY')
      print*,'year ',physPrefixForm,ichar
   end do
   
   ichar = INDEX(timestamp,'YYYY')
   print*,timestamp,ichar
   do while (ichar /= 0)
      write(timestamp(ichar:ichar+3),'(i4)') currYear
      ichar = INDEX(timestamp,'YYYY')
      print*,'year ',timestamp,ichar
   end do
   
   ichar = INDEX(physPrefixForm,'MM')
   do while (ichar /= 0)
      write(physPrefixForm(ichar:ichar+1),'(i2.2)') currMon
      ichar = INDEX(physPrefixForm,'MM')
   end do
      
   ichar = INDEX(physPrefixForm,'DD')
   do while (ichar /= 0)
      write(physPrefixForm(ichar:ichar+1),'(i2.2)') currDay   
      ichar = INDEX(physPrefixForm,'DD')
   end do 
   
   ichar = INDEX(physPrefixForm,'TSTSTSTS')
   do while (ichar /= 0)
      physPrefixForm = trim(physPrefixForm(:ichar-1))//trim(timestamp)//trim(physPrefixForm(ichar+8:))
      ichar = INDEX(physPrefixForm,'TSTSTSTS')
   end do
   
   ichar = INDEX(physPrefixForm,'GRIDX')
   tFile = trim(inDataDir)//trim(physPrefixForm(:ichar-1))//trim(tGridName)//trim(physPrefixForm(ichar+5:))//trim(fileSuffix)
   uFile = trim(inDataDir)//trim(physPrefixForm(:ichar-1))//trim(tGridName)//trim(physPrefixForm(ichar+5:))//trim(fileSuffix)
   vFile = trim(inDataDir)//trim(physPrefixForm(:ichar-1))//trim(tGridName)//trim(physPrefixForm(ichar+5:))//trim(fileSuffix)
   
   print*,tFile
    
   !dataprefix='xxxx/ORCA0083-N06_xxxxxxxx'
   !write(dataprefix(1:4),'(i4)')   currYear
   !write(dataprefix(19:26),'(i4,i2.2,i2.2)') currYear,currMon,currDay
   !fieldFile = trim(inDataDir)//'means/'//trim(dataprefix)//'d05'
   !medfieldFile = trim(inDataDir)//'medusa/'//trim(dataprefix)//'d05'
   !fieldFile = trim(inDataDir)//'means/2000/ORCA0083-N01_20000105d05'
   
   
   ! Read SSH
   if (readSSH .or. vvl) then
      hs(:,     :, nsp) = get2DfieldNC(trim(tFile), ssh_name)
      hs(imt+1, :, nsp) = hs(1,:,nsp)
   end if
   
   ! Depth at U, V, T points as 2D arrays
   allocate ( abyst(imt, jmt) , abysu(imt, jmt) , abysv(imt, jmt) )
   
   abyst = sum(dzt0(:,:,:), dim=3)
   abysu = sum(dzu(:,:,:,1), dim=3)
   abysv = sum(dzv(:,:,:,1), dim=3)
   
   if (vvl) then
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
   end if
   
   ! Read temperature 
   if (readTS) then
      xxx(:,:,:) = get3DfieldNC(trim(tFile), temp_name)
      tem(:,:,:,nsp) = xxx(:,:,km:1:-1)
      
      ! Read salinity
      xxx(:,:,:) = get3DfieldNC(trim(tFile), salt_name)
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
   end if     
   
   !
   ! Read u, v
   !
   uvel = get3DfieldNC(trim(uFile), trim(ueul_name))
   vvel = get3DfieldNC(trim(vFile), trim(veul_name))
   
   if (readBio) then
      ! Put oxygen in salinity field
      xxx(:,:,:) = get3DfieldNC(trim(medfieldFile)//'D.nc', 'TPP3')
      sal(:,:,:,nsp) = xxx(:,:,km:1:-1)
      
      xxx(:,:,:) = get3DfieldNC(trim(medfieldFile)//'P.nc', 'DIN')
      rho(:,:,:,nsp) = xxx(:,:,km:1:-1)
   end if 
   
   
   !
   ! Calculate zonal and meridional volume flux
   !
   ! Weight by (1 + ssh / depth)
   ! This is only an approximation of what NEMO really does
   ! but is accurate within 1% 
   !
   
   if (vvl) then
      do k = 1, km
      do j = 1, jmt
      do i = 1, imt
         dzt(i,j,k,nsp) = dzt0(i,j,k) * zstot(i,j)
      end do
      end do
      end do
   end if
   
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



