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
  real(wp), dimension(:), allocatable :: hsp
  real(wp), dimension(:), allocatable :: r_diag
  real(wp), dimension(:,:), allocatable :: r_mat
  real(wp), dimension(:,:), allocatable :: h_mat
  real(wp), dimension(:,:), allocatable :: b_mat
  real(wp), dimension(:), allocatable :: xa, x_post
  real(wp), dimension(:,:), allocatable :: gain_k, s_mat
  real(wp), dimension(:,:), allocatable :: sf_map
  real(wp), dimension(:,:,:), allocatable :: post_flux
  integer :: n_grid, n_obs, n_time

  character(len=256), dimension(:), allocatable :: receptor_path
  real(wp), dimension(:), allocatable :: receptor_gas, receptor_bg
  integer, dimension(:), allocatable :: receptor_year, receptor_month, receptor_day, &
                                        receptor_hour, receptor_minute, receptor_second


  type(datetime) :: a

  ! --- NetCDF I/O Test Variables (MOVED HERE) ---
  real(wp) :: matrix_out(3,3), matrix_in(3,3)

  ! 1 --- Namelist declarations ---
  namelist /input_config/ receptor_filename, prior_path, prior_name, prior_gas, &
  & foot_lat, foot_lon, foot_name
  namelist /output_config/ hsp_output_file, sf_filename, post_flux_filename
  namelist /model_config/ add_bg, obs_err_sd, model_err_sd, prior_err_sd
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

  close(10)
  print *, "---------------------------------------------------"
  
  ! 2 --- Receptors ---
  call example_csv_reader(receptor_filename, &
                          receptor_year, receptor_month, receptor_day, &
                          receptor_hour, receptor_minute, receptor_second, &
                          receptor_path, receptor_gas, receptor_bg)       


  print *, "---------------------------------------------------"
  
  ! 3 --- Footprints ---
  print *, "Reading prior"
  ! TODO: Make a function to read the prior
  ! based on the time in the footprint.
  ! In this case, I'm using annual emissions
  ! But in case of monthly, daily or hourly
  ! emissions, it may be different
  call read_3d_netcdf(prior_path, prior_in, prior_name)
  
  print *, "Prior shape: ", shape(prior_in)
  print *, "Prior sum: ", sum(prior_in)

  print *, "---------------------------------------------------"

  ! 4 Footprints
  print *, "Reading footprints"

  if (allocated(hsp)) deallocate(hsp)

  n_obs = size(receptor_path)
  allocate(hsp(n_obs))
  hsp = 0.0_wp

  n_grid = size(prior_in, 1) * size(prior_in, 2)
  print *, "State vector size (N_grid): ", n_grid
  allocate(h_mat(n_obs, n_grid))
  h_mat = 0.0_wp

    do i = 1, size(receptor_path)

       print *, trim(receptor_path(i))

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

      hsp(i) = sum(foot_in * prior_in)
      h_mat(i, :) = reshape(sum(foot_in * prior_in, dim=3), [n_grid])

    end do

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
  print *, "R diagonal (variance): ", r_diag(1)
  print *, "R shape: ", shape(r_mat)

  ! 5c Kalman Inversion
  print *, "---------------------------------------------------"
  print *, "Performing Kalman Inversion..."
  
  ! Prior state vector (scaling factors = 1.0)
  allocate(xa(n_grid))
  xa = 1.0_wp
  
  ! Prior covariance matrix B (diagonal)
  allocate(b_mat(n_grid, n_grid))
  b_mat = 0.0_wp
  do i = 1, n_grid
    b_mat(i, i) = prior_err_sd**2
  end do
  
  ! S = H B H' + R
  ! gain_k = B H' inv(S)
  allocate(s_mat(n_obs, n_obs))
  s_mat = matmul(h_mat, matmul(b_mat, transpose(h_mat))) + r_mat
  call invert_matrix(s_mat)
  
  allocate(gain_k(n_grid, n_obs))
  gain_k = matmul(b_mat, matmul(transpose(h_mat), s_mat))
  
  ! x_post = xa + K * (y - H*xa)
  ! Note: hsp calc above is essentially H*xa when xa=1
  allocate(x_post(n_grid))
  x_post = xa + matmul(gain_k, (receptor_gas - hsp))
  
  print *, "Posterior scaling factors range: ", minval(x_post), " to ", maxval(x_post)
  print *, "Average scaling factor: ", sum(x_post) / n_grid

  ! 5d Save Posterior Outputs
  print *, "---------------------------------------------------"
  print *, "Saving posterior outputs..."
  
  ! Reshape scaling factors to 2D
  allocate(sf_map(size(prior_in, 1), size(prior_in, 2)))
  sf_map = reshape(x_post, shape(sf_map))
  
  print *, "Writing scaling factors to: " // trim(sf_filename)
  call write_2d_netcdf(sf_filename, sf_map, "scaling_factors")
  
  ! Calculate 3D posterior flux
  n_time = size(prior_in, 3)
  allocate(post_flux(size(prior_in, 1), size(prior_in, 2), n_time))
  
  do i = 1, n_time
    post_flux(:, :, i) = prior_in(:, :, i) * sf_map
  end do
  
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

! 7 netcdf
  print *, "---------------------------------------------------"
  if (run_netcdf_test) then
    print '("--- NetCDF I/O Test ---")'
    
    ! Generate a sample matrix and invert it
    matrix_out = reshape([ 2.0_wp, 1.0_wp, 1.0_wp, &
                           4.0_wp, -6.0_wp, 0.0_wp, &
                           -2.0_wp, 7.0_wp, 2.0_wp ], [3,3])
    call invert_matrix(matrix_out)

    ! --- NetCDF I/O ---
    print '("Writing inverted matrix to file: ", A)', trim(netcdf_filename)
    call write_2d_netcdf(netcdf_filename, matrix_out, "inverted_A")

    ! Initialize input matrix to zeros
    matrix_in = 0.0_wp
    print '("Reading matrix from file...")'
    call read_2d_netcdf(netcdf_filename, matrix_in, "inverted_A")

    ! Verify the result
    if (all(abs(matrix_out - matrix_in) < 1.0e-12_wp)) then
      print '("SUCCESS: Data read back matches data written.")'
    else
      print '("FAILURE: Data mismatch.")'
    end if
  else
    print *, "Skipping NetCDF I/O Test."
  end if

end program main
