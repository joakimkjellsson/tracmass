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
   !    ===  Set up the grid for ORCA0083 configuration ===
   ! =============================================================
   ! Subroutine for defining the grid of the ORCA0083 config. 
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
   ! --- Where are coordinates stored, and horizontal grid, and vertical grid 
   !
   !coordFile = trim(topoDataDir)//'/2_INALT60.L120-KRS0020_mesh_mask.*'
   !hgridFile = trim(topoDataDir)//'/2_INALT60.L120-KRS0020_mesh_mask.*'
   !zgridFile = trim(topoDataDir)//'/2_INALT60.L120-KRS0020_mesh_mask.*'
   !bathyFile = trim(topoDataDir)//'/2_INALT60.L120-KRS0020_mesh_mask.*'
   !dx_name = 'e1t'
   !dy_name = 'e2t'
   !dxv_name = 'e1v'
   !dyu_name = 'e2u'
   !dz_1D_name = 'e3t_1d'
   !dzt_3D_name = 'e3t_0'
   !dzu_3D_name = 'e3u_0'
   !dzv_3D_name = 'e3v_0'
   !kBathy_name = 'mbathy'
   !gridIsUpsideDown = .true.
   !read3Ddz = .true. 
      
   !
   ! --- Read dx, dy at T points --- 
   !
   allocate ( e1t(imt,jmt) , e2t(imt,jmt) )
   e1t  = get2DfieldNC(trim(topoDataDir)//trim(hgridFile), dx_name)
   e2t  = get2DfieldNC(trim(topoDataDir)//trim(hgridFile), dy_name)
   dxdy(1:imt,1:jmt) = e1t(1:imt,1:jmt) * e2t(1:imt,1:jmt)
   deallocate ( e1t, e2t )
  
   !
   ! --- Read dy at U points and dx at V points --- 
   !
   dyu  = get2DfieldNC(trim(topoDataDir)//trim(hgridFile), dyu_name)
   dxv  = get2DfieldNC(trim(topoDataDir)//trim(hgridFile), dxv_name)
   dx   = dxv(imt/2, jmt/2)
   dy   = dyu(imt/2, jmt/2)
   
   !
   ! Read dz at T points without considering 
   ! bottom partial cells and variable volume  
   !
   dz = get1DfieldNC(trim(topoDataDir)//trim(zgridFile), dz_1D_name)
   if (gridIsUpsideDown) then
      do k=1,km
         kk=km+1-k
         dz(kk)=zlev(k)
         zlev(k)=zlev(k)+zlev(k-1)
      end do
   else
      do k=1,km
         dz(k)=zlev(k)
         zlev(k)=zlev(k)+zlev(k-1)
      end do
   end if
   
   !
   ! Read number of valid levels at U, V, T points
   ! as 2D array
   !
   kmt = get2DfieldNC(trim(topoDataDir)//trim(bathyFile), kBathy_name)
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
   
   ! Land-sea mask at north fold 
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
   allocate ( dzu(imt,jmt,km,2),dzv(imt,jmt,km,2), dzt0(imt,jmt,km) )
   
   if (read3Ddz) then
      
      print*,'  Reading 3D dz for u,v,t points ' 
      if (gridIsUpsideDown) then
         dzt0(:,:,km:1:-1)  = get3DfieldNC(trim(topoDataDir)//trim(zgridFile), dzt_3D_name)
         dzu(:,:,km:1:-1,1) = get3DfieldNC(trim(topoDataDir)//trim(zgridFile), dzu_3D_name)
         dzv(:,:,km:1:-1,1) = get3DfieldNC(trim(topoDataDir)//trim(zgridFile), dzv_3D_name)
      else
         dzt0(:,:,km:1:-1)  = get3DfieldNC(trim(topoDataDir)//trim(zgridFile), dzt_3D_name)
         dzu(:,:,km:1:-1,1) = get3DfieldNC(trim(topoDataDir)//trim(zgridFile), dzu_3D_name)
         dzv(:,:,km:1:-1,1) = get3DfieldNC(trim(topoDataDir)//trim(zgridFile), dzv_3D_name)
      end if
      
   else
      
      print*,'  Set dz horizontally constant for u,v,t points ' 
      print*,'  i.e. no partial steps '
      do j=1,jmt
         do i=1,imt
            dzt0(i,j,1:km) = dz(1:km)
            dzu(i,j,1:km,1)   = dz(1:km)
            dzv(i,j,1:km,1)   = dz(1:km)
         end do
      end do
      
   end if
      
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
      
   
end SUBROUTINE setupgrid
