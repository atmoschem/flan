program main
  use iso_fortran_env, only: wp => real64
  use linear_algebra, only: invert_matrix
  use io_manager, only: write_matrix_netcdf, read_matrix_netcdf
  ! Import the datetime module
  use datetime_module, only: datetime, datetime_type, now, strftime
  implicit none

  real(wp) :: matrix_out(3,3), matrix_in(3,3)
  character(len=20) :: filename = "matrix_data.nc"
  
  ! Variables for datetime test
  type(datetime_type) :: current_time
  character(len=100) :: formatted_time

  ! --- Datetime Test ---
  print '("--- Datetime Test ---")'
  call now(current_time)
  ! Format: YYYY-MM-DD HH:MM:SS
  call strftime(current_time, formatted_time, "%Y-%m-%d %H:%M:%S")
  print '("Current time is: ", A)', trim(formatted_time)

  !
  ! Generate a sample matrix and invert it
  matrix_out = reshape([ 2.0_wp, 1.0_wp, 1.0_wp, &
                         4.0_wp, -6.0_wp, 0.0_wp, &
                         -2.0_wp, 7.0_wp, 2.0_wp ], [3,3])
  call invert_matrix(matrix_out)

  !
  ! --- NetCDF Test ---
  print '("--- NetCDF I/O Test ---")'
  print '("Writing inverted matrix to file: ", A)', trim(filename)
  call write_matrix_netcdf(filename, matrix_out, "inverted_A")

  !
  ! Initialize input matrix to zeros to ensure it's populated by the read
  matrix_in = 0.0_wp
  print '("Reading matrix from file...")'
  call read_matrix_netcdf(filename, matrix_in, "inverted_A")

  !
  ! Verify the result
  if (all(abs(matrix_out - matrix_in) < 1.0e-12_wp)) then
    print '("SUCCESS: Data read back matches data written.")'
  else
    print '("FAILURE: Data mismatch.")'
  end if

end program main