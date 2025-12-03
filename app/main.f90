program main
  use iso_fortran_env, only: wp => real64
  use linear_algebra, only: invert_matrix
  use io_manager, only: write_matrix_netcdf, read_matrix_netcdf, example_csv_reader
  use datetime_module, only: datetime
  implicit none

  integer :: ios ! For I/O status

! --- Namelist Variables (with defaults) ---
  logical :: run_datetime_test = .true.
  logical :: run_netcdf_test   = .true.
  character(len=100) :: netcdf_filename = "matrix_data.nc"
  character(len=100) :: receptor_filename = "receptors.csv"


  character(len=100), dimension(:), allocatable :: receptor_path
  real(wp), dimension(:), allocatable :: receptor_gas, receptor_bg
  integer, dimension(:), allocatable :: receptor_year, receptor_month, receptor_day, &
                                        receptor_hour, receptor_minute, receptor_second


  type(datetime) :: a

  ! --- NetCDF I/O Test Variables (MOVED HERE) ---
  real(wp) :: matrix_out(3,3), matrix_in(3,3)

  ! 1 --- Namelist declarations ---
  namelist /receptor_config/ receptor_filename
  namelist /flan_config/ run_datetime_test, run_netcdf_test, netcdf_filename


! --- START OF EXECUTABLE CODE ---
  print *, "Reading configuration from 'namelists/config.nml'..."
  
  open(unit=10, file="namelists/config.nml", status="old", action="read", iostat=ios)

  if (ios /= 0) then
    print *, "Warning: Could not open 'namelists/config.nml'. Using default values."
  else
    ! Read the namelist group
    read(unit=10, nml=receptor_config, iostat=ios)
    if (ios /= 0) then
      print *, "Warning: Error reading namelist. Using defaults."
    end if
    close(10)
  end if

  ! 2 --- Receptors ---
  ! --- The bad line "a = datetime(...)" has been REMOVED from here ---
  call example_csv_reader(receptor_filename, &
                          receptor_year, receptor_month, receptor_day, &
                          receptor_hour, receptor_minute, receptor_second, &
                          receptor_path, receptor_gas, receptor_bg)       

  print *, "Configuration loaded."


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
    call write_matrix_netcdf(trim(netcdf_filename), matrix_out, "inverted_A")

    ! Initialize input matrix to zeros
    matrix_in = 0.0_wp
    print '("Reading matrix from file...")'
    call read_matrix_netcdf(trim(netcdf_filename), matrix_in, "inverted_A")

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