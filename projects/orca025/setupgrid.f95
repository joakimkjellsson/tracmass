SUBROUTINE setupgrid
  
   USE mod_precdef
   USE netcdf
   USE mod_param
   USE mod_vel
   
   USE mod_time
   USE mod_grid
   USE mod_name
   USE mod_vel
   USE mod_getfile
   
   IMPLICIT none
   ! =============================================================
   !    ===  Set up the grid for ORCA025 configuration ===
   ! =============================================================
   ! Subroutine for defining the grid of the ORCA025 config. 
   ! Run once before the loop starts.
   ! -------------------------------------------------------------
   ! The following arrays will be populated:
   !
   !  dxdy - Horizontal area of cells (T points)
   !  dz   - Thickness of standard level (T point) 
   !  dzt  - Time-invariant thickness of level (T point)
   !  dzu  - Time-invariant thickness of level (U point)
   !  dzv  - Time-invariant thickness of level (V point)
   !  kmt  - Number of levels from surface to seafloor (T point)
   !  kmu  - Number of levels from surface to seafloor (U point)
   !  kmv  - Number of levels from surface to seafloor (V point)
   !
   ! -------------------------------------------------------------
    
   ! === Init local variables for the subroutine ===
   INTEGER                                      :: i ,j ,k, n, kk, ii, &
   &                                               ip, jp, im, jm !! Loop indices
   REAL(DP), SAVE, ALLOCATABLE, DIMENSION(:,:)  :: e1t,e2t        !! dx, dy [m]
   REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:)    :: tmp4D
   CHARACTER (len=200)                          :: gridFile 
   
   
   map2D    = [3, 4,  1, 1 ]
   map3D    = [2, 3,  4, 1 ]
   ncTpos   = 1
   
   !
   ! --- Read dx, dy at T points --- 
   !
   print*,'setupgrid: set up dx,dy ',imt,jmt
   allocate ( e1t(imt+2,jmt) , e2t(imt+2,jmt) )
   gridFile = trim(inDataDir)//'domain/mesh_hgr.nc'
   print*,'read e1t, e2t from '//trim(gridFile)
   print*,'shape of e1t',shape(e1t),shape(e2t)
   e1t  = get2DfieldNC(trim(gridFile), 'e1t')
   print*,'e1t',e1t,size(e1t,1),size(e1t,2)
   e2t  = get2DfieldNC(trim(gridFile), 'e2t')
   print*,'e2t',e2t,size(e2t,1),size(e2t,2)
   dxdy(1:imt,1:jmt) = e1t(1:imt,1:jmt) * e2t(1:imt,1:jmt)
   deallocate ( e1t, e2t )
  
   !
   ! --- Read dy at U points and dx at V points --- 
   !
   dyu  = get2DfieldNC(gridFile, 'e2u')
   dxv  = get2DfieldNC(gridFile, 'e1v')
   dx   = dxv(imt/2, jmt/2)
   dy   = dyu(imt/2, jmt/2)
   
   !
   ! Read dz at T points without considering 
   ! bottom partial cells and variable volume  
   !
   print*,'setupgrid: set up dz'
   gridFile = trim(inDataDir)//'domain/mesh_zgr.nc'
   dz = get1DfieldNC(gridFile, 'e3t_0')
   dz(1:km) = dz(km:1:-1)
   !do k=1,km
   !   kk=km+1-k
   !   dz(kk)=zlev(k)
   !   zlev(k)=zlev(k)+zlev(k-1)
   !end do
   
   !
   ! Read number of valid levels at U, V, T points
   ! as 2D array
   !
   print*,'setupgrid: set up bathymetry'
   kmt = get2DfieldNC(gridFile, 'mbathy')
   allocate ( kmu(imt,jmt), kmv(imt,jmt) )
   
   kmu=0 ; kmv=0
   do j=1,jmt
      jp=j+1
      if(jp == jmt+1) jp=jmt
      do i=1,imt
         ip=i+1
         if(ip == imt+1) ip=1
         kmu(i,j)=min(kmt(i,j), kmt(ip,j),KM)
         kmv(i,j)=min(kmt(i,j), kmt(i,jp),KM)
      enddo
   enddo

   do i=4, imt
      ii = imt + 4 - i
      kmv(i,jmt) = kmv(ii,jmt-3)
   enddo
  
   !
   ! Read layer thickness at U, V, T points 
   ! without considering variable volume.
   !
   ! SSH variability is accounted for in readfield each time step
   !
   print*,'setupgrid: set up layer thickness'
   allocate ( dzu(imt,jmt,km,2),dzv(imt,jmt,km,2), dzt0(imt,jmt,km) )
   
   dzt0(:,:,:) = get3DfieldNC(gridFile, 'e3t')
   dzu(:,:,:,1) = get3DfieldNC(gridFile, 'e3u')
   dzv(:,:,:,1) = get3DfieldNC(gridFile, 'e3v')
   
   !
   ! Ensure thickness is zero in invalid points
   !
   do n=1,2
      do k=1,km
         where (k > kmt(1:imt,1:jmt))
            dzt0(:,:,k) = 0
         end where
         where (k > kmu(1:imt,1:jmt))
            dzu(:,:,k,n) = 0
         end where
         where (k > kmv(1:imt,1:jmt))
            dzv(:,:,k,n) = 0
         end where
      enddo 
   end do
   
   ! Reverse grid 
   allocate( tmp4D(imt,jmt,km,2) )
   
   tmp4D(1:imt,1:jmt,1:km,1) = dzt0(1:imt,1:jmt,1:km)
   do k=1,km
      dzt0(:,:,k) = tmp4D(:,:,km+1-k,1) 
   end do
   
   tmp4D(1:imt,1:jmt,1:km,:) = dzu(1:imt,1:jmt,1:km,:)
   do k=1,km
      dzu(:,:,k,:) = tmp4D(:,:,km+1-k,:) 
   end do
   
   tmp4D(1:imt,1:jmt,1:km,:) = dzv(1:imt,1:jmt,1:km,:)
   do k=1,km
      dzv(:,:,k,:) = tmp4D(:,:,km+1-k,:) 
   end do
   
   deallocate( tmp4d )
   
end SUBROUTINE setupgrid
