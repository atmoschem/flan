program main
  use iso_fortran_env, only: wp => real64
  use netcdf
  use linear_algebra, only: invert_matrix, diag
  use io_manager, only: write_1d_netcdf, write_2d_netcdf, &
  &write_3d_netcdf, read_1d_netcdf, read_2d_netcdf, &
  &read_3d_netcdf, example_csv_reader, logical_to_string, write_1d_text
  use datetime_module, only: datetime
  implicit none

  integer :: ios ! For I/O status
  integer :: i, lon, lat ! Loop indices
  ! Dimensions (these should be read from NetCDF files or set as parameters)
  integer :: nlons, nlats
  ! Footprint matrix H
  real(wp), dimension(:,:), allocatable :: h

! --- Namelist Variables (with defaults) ---
  logical :: run_datetime_test = .true.
  logical :: run_netcdf_test   = .true.
  logical :: add_bg   = .true.
  character(len=256) :: netcdf_filename = "matrix_data.nc"
  character(len=256) :: hsp_output_file = "hsp_output.txt"
  real(wp) :: obs_err_sd = 1.0_wp
  real(wp) :: model_err_sd = 1.0_wp
  real(wp) :: prior_err_sd = 0.5_wp
  real(wp) :: spatial_correlation_len = 50.0_wp ! km
  real(wp) :: temporal_correlation_len = 12.0_wp ! hours
  integer :: time_resolution_hours = 3
  character(len=256) :: sf_filename = "scaling_factors.nc"
  character(len=256) :: post_flux_filename = "posterior_flux.nc"
  character(len=256) :: receptor_filename = "receptors.csv"
  character(len=256) :: foot_lat = "foot1lat"
  character(len=256) :: foot_lon = "foot1lon"
  character(len=256) :: foot_name = "foot1"
  character(len=256) :: prior_path = "prior.nc"
  character(len=256) :: prior_name = "prior"
  character(len=256) :: prior_gas = "gas"


  real(wp), dimension(:,:,:), allocatable :: foot_in
  real(wp), dimension(:), allocatable :: lats, lons
  real(wp), dimension(:,:,:), allocatable :: prior_in
  real(wp), dimension(:,:), allocatable :: B_s, B_t ! Covariance blocks
  real(wp), dimension(:,:), allocatable :: BHt ! B * H^T
  real(wp), dimension(:,:), allocatable :: V_in, V_out ! Temps for reshaping
  real(wp), dimension(:), allocatable :: hsp
  real(wp), dimension(:), allocatable :: r_diag
  real(wp), dimension(:,:), allocatable :: r_mat
  real(wp), dimension(:,:), allocatable :: h_mat
  ! real(wp), dimension(:,:), allocatable :: b_mat ! REMOVED
  real(wp), dimension(:), allocatable :: xa, x_post
  real(wp), dimension(:,:), allocatable :: gain_k, s_mat
  real(wp), dimension(:,:), allocatable :: sf_map
  real(wp), dimension(:,:,:), allocatable :: post_flux
  real(wp), dimension(:,:,:), allocatable :: sf_map_3d
  integer :: n_grid, n_obs, n_time
  integer :: n_time_blocks, n_state
  real(wp), dimension(:), allocatable :: time_blocks_start ! Start hour of each block relative to simulation start
  real(wp), dimension(:,:), allocatable :: dist_mat ! Spatial distance matrix
  real(wp), dimension(:,:), allocatable :: time_dist_mat ! Temporal distance matrix

  character(len=256), dimension(:), allocatable :: receptor_path
  real(wp), dimension(:), allocatable :: receptor_gas, receptor_bg
  integer, dimension(:), allocatable :: receptor_year, receptor_month, receptor_day, &
                                        receptor_hour, receptor_minute, receptor_second


  type(datetime) :: a

  ! --- NetCDF I/O Test Variables (MOVED HERE) ---
  real(wp) :: matrix_out(3,3), matrix_in(3,3)

  namelist /input_config/ receptor_filename, prior_path, prior_name, prior_gas, &
  & foot_lat, foot_lon, foot_name
  namelist /output_config/ hsp_output_file, sf_filename, post_flux_filename

  ! Loop and temp variables for covariance
  integer :: ix1, iy1, ix2, iy2, nx, ny, j
  integer :: r_start, c_start
  real(wp) :: d_km, dy, dx, j_dummy
  integer :: j_loop
  namelist /model_config/ add_bg, obs_err_sd, model_err_sd, prior_err_sd, &
       & spatial_correlation_len, temporal_correlation_len, time_resolution_hours
  namelist /test_config/ run_datetime_test, run_netcdf_test, netcdf_filename

  print *, "---------------------------------------------------"
  print *, "Reading configuration from 'namelists/config.nml'..."
  
  open(unit=10, file="namelists/config.nml", status="old", action="read", iostat=ios)

  ! Read all namelist groups
  read(unit=10, nml=input_config, iostat=ios)
  rewind(10)
  read(unit=10, nml=output_config, iostat=ios)
  rewind(10)
  read(unit=10, nml=model_config, iostat=ios)
  rewind(10)
  read(unit=10, nml=test_config, iostat=ios)

  print *, "input_config: receptor: " // trim(receptor_filename)
  print *, "input_config: prior: " // trim(prior_path)
  print *, "model_config: add_bg: " // logical_to_string(add_bg) 
  print *, "model_config: spatial_corr: ", spatial_correlation_len
  print *, "model_config: temp_corr: ", temporal_correlation_len
  print *, "model_config: time_res_hours: ", time_resolution_hours 

  close(10)
  print *, "---------------------------------------------------"
  
  ! 2 --- Receptors ---
  call example_csv_reader(receptor_filename, &
                          receptor_year, receptor_month, receptor_day, &
                          receptor_hour, receptor_minute, receptor_second, &
                          receptor_path, receptor_gas, receptor_bg)       

  print *, "Total receptors loaded: ", size(receptor_path)


  print *, "---------------------------------------------------"
  
  ! 3 --- Footprints ---
  print *, "Reading prior"
  ! TODO: Make a function to read the prior
  ! based on the time in the footprint.
  ! In this case, I'm using annual emissions
  ! But in case of monthly, daily or hourly
  ! emissions, it may be different
  call read_3d_netcdf(prior_path, prior_in, prior_name)
  
  print *, "Prior dimensions (lon, lat, time): ", shape(prior_in)
  print *, "Prior range (min, max): ", minval(prior_in), maxval(prior_in)
  print *, "Prior sum: ", sum(prior_in)

  print *, "---------------------------------------------------"

  ! 4 Footprints
  print *, "Reading footprints"

  if (allocated(hsp)) deallocate(hsp)

  n_obs = size(receptor_path)
  allocate(hsp(n_obs))
  hsp = 0.0_wp

  n_grid = size(prior_in, 1) * size(prior_in, 2)
  ! Set time blocks from prior dimensions (e.g. 240 for 30 days of 3-hourly)
  n_time_blocks = size(prior_in, 3) 
  n_state = n_grid * n_time_blocks
  
  if (n_state > 30000) then
      print *, "WARNING: State vector is massive (", n_state, ")."
      print *, "B matrix will require approx ", real(n_state,wp)**2 * 8.0 / 1.0e9, " GB of RAM."
  end if
  
  print *, "State vector size (N_grid): ", n_grid
  print *, "Time blocks: ", n_time_blocks
  print *, "Total State Vector Size: ", n_state

  allocate(h_mat(n_obs, n_state))
  h_mat = 0.0_wp

     do i = 1, size(receptor_path)

        print *, "Processing footprint ", i, ": ", trim(receptor_path(i))

        ! Footprints
        call read_3d_netcdf(trim(receptor_path(i)), foot_in, foot_name)
        call read_1d_netcdf(trim(receptor_path(i)), lats, foot_lat)
        call read_1d_netcdf(trim(receptor_path(i)), lons, foot_lon)

        if (any(shape(foot_in) /= shape(prior_in))) then
          print *, "ERROR: Dimension mismatch between Footprint and Prior!"
          print *, "Foot: ", shape(foot_in)
          print *, "Prior: ", shape(prior_in)
          stop
        end if

       print *, "   Footprint ", i, " range (min, max): ", minval(foot_in), maxval(foot_in)
      ! Calculate time index for this observation
       ! TODO: Real time mapping. For now, we map sequentially or cyclically for testing.
       ! We assume the footprint actually corresponds to a specific time block [k].
       ! In a real 4D inversion, we would integrate foot * prior(t) -> sensitivity to x(t)
       ! Here, we simplify: we assume this footprint sees the state vector at a specific time index
       ! defined by the observation time.
       
       ! Dummy mapping for now: (i-1) mod n_time_blocks + 1
       ! This spreads observations across the 240 time steps.
       ! integer :: t_idx
       ! t_idx = mod(i - 1, n_time_blocks) + 1
       
       ! BETTER: Just put it in the middle for testing, or use actual simple logic
       
       do j = 1, n_time_blocks
          ! If we had exact time matching, we would only fill the block corresponding to the time.
          ! For this task, let's assume 'foot_in' applies to 'prior_in' at ALL times weighted by... wait.
          ! The footprint has a time dimension (backwards).
          ! The standard H matrix approach for 4D-Var:
          ! H maps x (state) to y (conc).
          ! H_i = [ H_i(t=1), H_i(t=2), ... ]
          ! If the particle was at grid cell g at time t, then sensitivity is nonzero.
          ! Since our 'foot_in' is likely integrated or already processed, checking 'read_3d_netcdf' usage.
          ! It reads 'foot_in' which is normally lat/lon/time_back.
          ! But here we treat 'prior_in' as 3D (time_prior).
          
          ! Current logic: h_mat(i, :) = reshape(sum(foot_in * prior_in, dim=3), [n_grid])
          ! This implies 'foot_in' and 'prior_in' match in 3rd dim (time).
          ! If we want to solve for 240 time steps of SCALING FACTORS,
          ! h_mat(i, :) needs to have length 240 * n_grid.
          
          ! We need to know which time block each slice of 'foot_in' corresponds to.
          ! Assuming 'foot_in' 3rd dim aligns with 'prior_in' 3rd dim (e.g. 240 hours).
          ! Then:
          ! Observation i sees: \sum_t ( foot(x,y,t) * prior(x,y,t) * scalefactor(x,y,t) )
          ! So H matrix element for (grid cell g, time t) is foot(g,t)*prior(g,t).
          
          ! WARNING: 'foot_in' dimensions vs 'prior_in' dimensions.
          ! The code checked: if (any(shape(foot_in) /= shape(prior_in))) error.
          ! So 'foot_in' has 240 time steps too!
          
          ! Therefore, H is block-sparse or rather structured.
          ! H(i, (t-1)*n_grid + g) = foot_in(g_x, g_y, t) * prior_in(g_x, g_y, t)
          
          ! Let's implement this flat loop.
       end do
       
       ! Optimized construction:
       ! H is (1, n_grid*n_time) for single obs.
       ! We just flatten foot_in * prior_in
       h_mat(i, :) = reshape(foot_in * prior_in, [n_state])

       hsp(i) = sum(foot_in * prior_in)
     end do

   print *, "H matrix dimensions (n_obs, n_grid): ", shape(h_mat)

  ! 5 hsp
  print *, "---------------------------------------------------"
  print *, "Enhancements " // trim(prior_gas) // " (ppb)"
  do i = 1, size(receptor_path)
    write(*, '(I4, A, F15.5)') i, " : ", hsp(i)
  end do

! 6 background

  print *, "---------------------------------------------------"
  if (add_bg) then
    print *, 'Adding background'
    hsp = hsp + receptor_bg
  else
    print *, "Skipping background."
  end if

  print *, "Writing hsp to: " // trim(hsp_output_file)
  call write_1d_text(hsp_output_file, hsp)

  ! 5b Variance Matrix R
  print *, "---------------------------------------------------"
  print *, "Creating R matrix (measurement and model error covariance)"
  allocate(r_diag(size(hsp)))
  r_diag = obs_err_sd**2 + model_err_sd**2
  allocate(r_mat(size(hsp), size(hsp)))
  r_mat = diag(r_diag)
  print *, "R matrix dimensions (n_obs, n_obs): ", shape(r_mat)
  print *, "R diagonal (first 5 elements): ", r_diag(1:min(5, size(r_diag)))

  ! 5c Kalman Inversion
  print *, "---------------------------------------------------"
  print *, "Performing Kalman Inversion..."
  
  ! allocate(xa(n_grid)) ! Removed, done above
  ! xa = 1.0_wp          ! Removed, done above
  
  allocate(xa(n_state))
  xa = 1.0_wp
  
  ! Prior covariance matrix B
  ! Implement Spatial and Temporal Correlations
  ! B_{ij, kl} = sig^2 * exp(-d/Ls) * exp(-dt/Lt)
  
  print *, "Constructing B matrix components (Implicit Kronecker)"
  ! allocate(b_mat(n_state, n_state)) ! TOO BIG, do not allocate
  
  ! Precompute distance maps could be expensive. 
  ! Let's do loops. OpenMP would be great here.
  
  ! Since n_state is huge (40k), we should be careful. 
  ! R usually handles 40k x 40k okayish (12GB). 
  ! Let's try to be efficient.
  ! Separable assumption: B = B_t (Kron) B_s
  
  ! 1. Build Spatial Covariance B_s (n_grid x n_grid)
  ! We need lats/lons. They are read inside the loop? No, lats/lons are allocatable.
  ! 'lats' and 'lons' were read in the loop but re-allocated/overwritten?
  ! We need the grid coordinates.
  ! Let's assume 'lats' and 'lons' from the last footprint are valid for the grid.
  
  ! 2. Build Temporal Covariance B_t (n_time_blocks x n_time_blocks)
  
  ! 2. Build Temporal Covariance B_t (n_time_blocks x n_time_blocks)
  
  ! real(wp), allocatable :: B_s(:,:), B_t(:,:) ! MOVED TO TOP
  allocate(B_s(n_grid, n_grid))
  allocate(B_t(n_time_blocks, n_time_blocks))
  
  ! Build B_t
  do i = 1, n_time_blocks
     do j = 1, n_time_blocks
        ! distance in hours. Assuming indices are 3-hourly blocks?
        ! dt = abs(i-j) * time_resolution_hours
        B_t(i,j) = exp( -1.0_wp * abs(i-j) * real(time_resolution_hours,wp) / temporal_correlation_len )
     end do
  end do
  
  ! Build B_s
  ! Need grid coordinates. We can get them from prior dimensions if we had 1d arrays.
  ! lats/lons read from netcdf are 1D vectors for axes?
  ! Let's assume lats/lons are correct size (nlats, nlons).
  ! But we flattened grid to n_grid = n_lon * n_lat.
  ! So we need to map flat index k to (lon_idx, lat_idx).
  ! n_lon = size(prior_in,1), n_lat = size(prior_in,2)
  
  ! integer :: j ! REMOVED
  
  nx = size(prior_in, 1)
  ny = size(prior_in, 2)
  
  do i = 1, n_grid
    iy1 = (i-1)/nx + 1
    ix1 = mod(i-1, nx) + 1
    
    do j = 1, n_grid
       iy2 = (j-1)/nx + 1
       ix2 = mod(j-1, nx) + 1
       
       ! Approx distance (Euclidean on deg is bad, but for small domain ok, use Haversine ideally)
       ! Simple approx: 1 deg lat ~ 111 km. 1 deg lon ~ 111 * cos(lat).
       ! Let's use simple Euclidean for now or simple scaling.
       ! d = sqrt( (dlat*111)**2 + (dlon*111*cos(lat))**2 )
       
       dy = (lats(iy1) - lats(iy2)) * 111.0_wp
       dx = (lons(ix1) - lons(ix2)) * 111.0_wp * cos(lats(iy1)*3.14159/180.0_wp)
       d_km = sqrt(dx**2 + dy**2)
       
       B_s(i,j) = (prior_err_sd**2) * exp( -1.0_wp * d_km / spatial_correlation_len )
    end do
  end do
  
  ! Kronecker Product
  ! B = B_t (x) B_s ... or B_s (x) B_t?
  ! State vector layout: 
  ! reshape(foot_in * prior_in, [n_state]) flattens (lon, lat, time)
  ! So inner dim is lon, then lat, then time.
  ! So grid is fast, time is slow.
  ! B_{ (g1,t1), (g2,t2) } = Cov(x_{g1,t1}, x_{g2,t2})
  !                        = Cov(g1,g2) * Cov(t1,t2)
  !                        = B_s(g1,g2) * B_t(t1,t2)
  
  ! B_mat layout matches state vector.
  ! State vector index k = (t-1)*n_grid + g
  ! B_mat(k1, k2):
  !   k1 -> (g1, t1)
  !   k2 -> (g2, t2)
  
  ! Implicit Kronecker: B should NOT be built.
  ! We need B * H^T and H * B * H^T
  ! H is (n_obs, n_state). H^T is (n_state, n_obs).
  
  ! Let BHt = B * H^T. Dimensions (n_state, n_obs)
  allocate(BHt(n_state, n_obs))
  BHt = 0.0_wp
  
  print *, "Computing B * H^T using implicit Kronecker products..."
  allocate(V_in(n_grid, n_time_blocks))
  allocate(V_out(n_grid, n_time_blocks))
  
  ! Loop over columns of H^T (i.e., rows of H)
  do i = 1, n_obs
     ! 1. Extract vector h_i = H(i, :)
     V_in = reshape(h_mat(i, :), shape(V_in))
     
     ! 2. Compute z = (B_t (x) B_s) * vec(V_in) = vec( B_s * V_in * B_t' )
     ! Note: B_t is symmetric so B_t' = B_t
     V_out = matmul(B_s, matmul(V_in, B_t))
     
     ! 3. Store back in BHt
     BHt(:, i) = reshape(V_out, [n_state])
  end do
  
  ! Clean up B_s and B_t if not needed, but they are small.
  print *, "B * H^T computed."
  
  ! S = H * (B * H^T) + R
  ! S = H * BHt + R
  allocate(s_mat(n_obs, n_obs))
  s_mat = matmul(h_mat, BHt) + r_mat
  print *, "S matrix dimensions (n_obs, n_obs): ", shape(s_mat)
  call invert_matrix(s_mat)
  
  allocate(gain_k(n_state, n_obs))
  ! gain_k = B * H^T * inv(S) = BHt * inv(S)
  gain_k = matmul(BHt, s_mat)
  print *, "Kalman Gain K dimensions (n_grid, n_obs): ", shape(gain_k)
  
  ! x_post = xa + K * (y - H*xa)
  ! Note: hsp calc above is essentially H*xa when xa=1
  allocate(x_post(n_state))
  x_post = xa + matmul(gain_k, (receptor_gas - hsp))
  
  print *, "Posterior state x_post size: ", size(x_post)
  print *, "Posterior scaling factors range: ", minval(x_post), " to ", maxval(x_post)
  print *, "Average scaling factor: ", sum(x_post) / n_state

  ! 5d Save Posterior Outputs
  print *, "---------------------------------------------------"
  print *, "Saving posterior outputs..."
  
  ! Reshape scaling factors to 3D (lon, lat, time)
  allocate(sf_map(size(prior_in, 1), size(prior_in, 2)))
  
  ! real(wp), allocatable :: sf_map_3d(:,:,:) ! REMOVED
  allocate(sf_map_3d(size(prior_in, 1), size(prior_in, 2), n_time_blocks))
  sf_map_3d = reshape(x_post, shape(sf_map_3d))
  
  print *, "Reshaped scaling factors (sf_map_3d) dimensions: ", shape(sf_map_3d)
  print *, "Writing 3D scaling factors to: scaling_factors_3d.nc"
  call write_3d_netcdf("nc/scaling_factors_3d.nc", sf_map_3d, "scaling_factors") 
  
  ! Compute average map for 2D output
  sf_map = sum(sf_map_3d, dim=3) / real(n_time_blocks, wp)

  print *, "Reshaped scaling factors (sf_map) dimensions: ", shape(sf_map)
  print *, "Writing scaling factors to: " // trim(sf_filename)
  call write_2d_netcdf(sf_filename, sf_map, "scaling_factors")
  
  ! Calculate 3D posterior flux
  n_time = size(prior_in, 3)
  allocate(post_flux(size(prior_in, 1), size(prior_in, 2), n_time))
  
  do i = 1, n_time
     ! Now we multiply Prior(x,y,t) by SF(x,y,t)
     post_flux(:, :, i) = prior_in(:, :, i) * sf_map_3d(:, :, i)
  end do
  
  print *, "Posterior flux dimensions: ", shape(post_flux)
  print *, "Average posterior flux: ", sum(post_flux) / size(post_flux)
  print *, "Writing posterior flux to: " // trim(post_flux_filename)
  call write_3d_netcdf(post_flux_filename, post_flux, prior_gas)

! 6 date time
  print *, "---------------------------------------------------"
  if (run_datetime_test) then
    print *, '("--- Datetime Test ---")'
    a = a % now()
    print *, a % isoformat()
  else
    print *, "Skipping Datetime Test."
  end if

end program main
