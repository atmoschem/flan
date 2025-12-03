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

  real(wp), dimension(:,:,:), allocatable :: foot_in
  real(wp), dimension(:), allocatable :: lats, lons

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
  namelist /prior_config/ prior_path


! --- START OF EXECUTABLE CODE ---
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

  ! 2 --- Receptors ---
  call example_csv_reader(receptor_filename, &
                          receptor_year, receptor_month, receptor_day, &
                          receptor_hour, receptor_minute, receptor_second, &
                          receptor_path, receptor_gas, receptor_bg)       


  ! 4 Footprints
    print *, "Reading footprints"

    do i = 1, size(receptor_path)

       print *, trim(receptor_path(i))

       ! Footprints
       call read_3d_netcdf(trim(receptor_path(i)), foot_in, foot_name)
       call read_1d_netcdf(trim(receptor_path(i)), lats, foot_lat)
       call read_1d_netcdf(trim(receptor_path(i)), lons, foot_lon)

!       h(i, :) =  reshape(foot_in, [size(lons)* size(lats)])

    end do
    ! TODO: Add integer-to-string conversion function
    ! print *, "h (", size(h, 1), ",", size(h, 2), ")"




  if (run_datetime_test) then
    print *, '("--- Datetime Test ---")'
    a = a % now()
    print *, a % isoformat()
  else
    print *, "Skipping Datetime Test."
  end if

  ! --- 2. Conditional NetCDF I/O Test ---
 
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