program main
  use iso_fortran_env, only: wp => real64
  use netcdf
  use linear_algebra, only: invert_matrix
  use io_manager, only: write_1d_netcdf, write_2d_netcdf, &
  &write_3d_netcdf, read_1d_netcdf, read_2d_netcdf, &
  &read_3d_netcdf, example_csv_reader
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
  character(len=256) :: netcdf_filename = "matrix_data.nc"
  character(len=256) :: receptor_filename = "receptors.csv"
  character(len=256) :: foot_lat = "foot1lat"
  character(len=256) :: foot_lon = "foot1lon"
  character(len=256) :: foot_name = "foot1"
  character(len=256) :: prior_path = "prior.nc"
  character(len=256) :: prior_name = "prior"


  real(wp), dimension(:,:,:), allocatable :: foot_in
  real(wp), dimension(:), allocatable :: lats, lons
  real(wp), dimension(:,:,:), allocatable :: prior_in
  real(wp), dimension(:), allocatable :: hsp

  character(len=256), dimension(:), allocatable :: receptor_path
  real(wp), dimension(:), allocatable :: receptor_gas, receptor_bg
  integer, dimension(:), allocatable :: receptor_year, receptor_month, receptor_day, &
                                        receptor_hour, receptor_minute, receptor_second


  type(datetime) :: a

  ! --- NetCDF I/O Test Variables (MOVED HERE) ---
  real(wp) :: matrix_out(3,3), matrix_in(3,3)

  ! 1 --- Namelist declarations ---
  namelist /flan_config/ run_datetime_test, run_netcdf_test, netcdf_filename
  namelist /receptor_config/ receptor_filename
  namelist /footprint_config/ foot_lat, foot_lon, foot_name
  namelist /prior_config/ prior_path, prior_name

  print *, "---------------------------------------------------"
  print *, "Reading configuration from 'namelists/config.nml'..."
  
  open(unit=10, file="namelists/config.nml", status="old", action="read", iostat=ios)

 ! Read the namelist group
  read(unit=10, nml=receptor_config, iostat=ios)
  print *, "Receptor filename: " // receptor_filename

  rewind(10)
  read(10, nml=footprint_config, iostat=ios)
  print *, "Footprint lat name: " // foot_lat
  print *, "Footprint lon name: " // foot_lon
  print *, "Footprint name: " // foot_name

  rewind(10)
  read(10, nml=prior_config, iostat=ios)
  print *, "Prior NetCDF filename: " // prior_path

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

  allocate(hsp(size(receptor_path)))

  hsp = 0.0_wp

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

    end do

  ! 5 hsp
  print *, "---------------------------------------------------"
  print *, "Enhancements"
  do i = 1, size(receptor_path)
      print *, i, ": ", hsp(i)
  end do

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