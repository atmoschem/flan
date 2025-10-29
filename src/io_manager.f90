module io_manager
  use iso_fortran_env, only: wp => real64
  use netcdf
  use csv_module
  use stdlib_error, only: error_stop
  implicit none

  ! CORRECT: Make both subroutines public
  public :: write_matrix_netcdf, read_matrix_netcdf
  public :: example_csv_reader

contains

  subroutine check(status)
   ! Helper subroutine for error checking
    integer, intent(in) :: status
    if (status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "NetCDF Error"
    end if
  end subroutine check

  subroutine write_matrix_netcdf(filename, data_matrix, varname)
    character(len=*), intent(in) :: filename
    real(wp), intent(in) :: data_matrix(:,:)
    character(len=*), intent(in) :: varname

    integer :: ncid, dim1_id, dim2_id, var_id
    integer :: status
    integer :: dim1_size, dim2_size

    dim1_size = size(data_matrix, 1)
    dim2_size = size(data_matrix, 2)

    status = nf90_create(filename, nf90_clobber, ncid)
    call check(status)

    status = nf90_def_dim(ncid, "dim1", dim1_size, dim1_id)
    call check(status)
    status = nf90_def_dim(ncid, "dim2", dim2_size, dim2_id)
    call check(status)

    status = nf90_def_var(ncid, varname, nf90_double, [dim2_id, dim1_id], var_id)
    call check(status)

    status = nf90_enddef(ncid)
    call check(status)

    status = nf90_put_var(ncid, var_id, data_matrix)
    call check(status)

    status = nf90_close(ncid)
    call check(status)
  end subroutine write_matrix_netcdf

 ! CORRECT: The missing subroutine is now included.
  subroutine read_matrix_netcdf(filename, data_matrix, varname)
    character(len=*), intent(in) :: filename
    real(wp), intent(inout) :: data_matrix(:,:)
    character(len=*), intent(in) :: varname

    integer :: ncid, var_id, status

   ! 1. Open the file for reading.
    status = nf90_open(filename, nf90_nowrite, ncid)
    call check(status)

   ! 2. Get the ID of the variable by its name.
    status = nf90_inq_varid(ncid, varname, var_id)
    call check(status)

   ! 3. Read the data from the variable into the array.
    status = nf90_get_var(ncid, var_id, data_matrix)
    call check(status)

   ! 4. Close the file.
    status = nf90_close(ncid)
    call check(status)
  end subroutine read_matrix_netcdf

! [This code replaces the existing 'example_csv_reader' subroutine in src/io_manager.f90]

! --- Example CSV Reader using csv-fortran ---
  subroutine example_csv_reader(filename)
    character(len=*), intent(in) :: filename
    type(csv_file) :: csv
    character(len=:), allocatable :: header(:), line_data(:)
    logical :: at_end
    integer :: i, ierr
    
    ! Variables to hold column indices
    integer :: time_idx = -1, gas_idx = -1, year_idx = -1, month_idx = -1, &
               day_idx = -1, hour_idx = -1, minute_idx = -1, path_idx = -1

    ! Variables to hold data from a row
    character(len=100) :: time_str, gas_str, path_str
    character(len=4)   :: year_str
    character(len=2)   :: month_str, day_str, hour_str, minute_str
    integer :: year_int, month_int, day_int, hour_int, minute_int
    character(len=100) :: output_filename
    
    ! Column names to find
    character(len=*), parameter :: TIME_COL = "time_col"
    character(len=*), parameter :: GAS_COL = "gas_col"
    character(len=*), parameter :: YEAR_COL = "year"
    character(len=*), parameter :: MONTH_COL = "month"
    character(len=*), parameter :: DAY_COL = "day"
    character(len=*), parameter :: HOUR_COL = "hour"
    character(len=*), parameter :: MINUTE_COL = "minute"
    character(len=*), parameter :: PATH_COL = "path_col"

    ! Helper function to find column index
    contains
      function find_col_index(header_array, col_name) result(idx)
        character(len=:), allocatable, intent(in) :: header_array(:)
        character(len=*), intent(in) :: col_name
        integer :: idx, j
        idx = -1 ! Default if not found
        do j = 1, size(header_array)
          if (trim(header_array(j)) == col_name) then
            idx = j
            exit
          end if
        end do
      end function find_col_index

    ! --- Main Subroutine Logic ---

    ! Open the CSV file
    call csv%open(filename, n_cols = 8, status_ok = at_end)
!    if (.not. csv%is_open()) then
!      call error_stop("Error: Could not open CSV file: " // trim(filename))
!    end if

    ! 1. Read header to find column indices
    call read_line(csv, header, at_end)
    if (at_end) then
      print *, "Warning: CSV file is empty."
      call close(csv)
      return
    end if
    
    ! Find indices for all required columns
    time_idx   = find_col_index(header, TIME_COL)
    gas_idx    = find_col_index(header, GAS_COL)
    year_idx   = find_col_index(header, YEAR_COL)
    month_idx  = find_col_index(header, MONTH_COL)
    day_idx    = find_col_index(header, DAY_COL)
    hour_idx   = find_col_index(header, HOUR_COL)
    minute_idx = find_col_index(header, MINUTE_COL)
    path_idx   = find_col_index(header, PATH_COL)
    
    ! Basic check if all columns were found
    if (minval([time_idx, gas_idx, year_idx, month_idx, day_idx, hour_idx, minute_idx, path_idx]) < 1) then
       print *, "Error: Could not find one or more required columns in CSV header."
       print *, "Found indices: ", time_idx, gas_idx, year_idx, month_idx, day_idx, hour_idx, minute_idx, path_idx
       call close_csv(csv)
       call error_stop("Missing CSV columns")
    end if
    
    print *, "Successfully mapped all columns."

    ! 2. Loop through data rows
    do
      call read_line(csv, line_data, at_end)
      if (at_end) exit ! Exit loop at end of file
      
      ! 3. Extract string data from the correct columns
      time_str = line_data(time_idx)
      gas_str  = line_data(gas_idx)
      path_str = line_data(path_idx)
      
      ! 4. Convert numeric strings to integers for processing
      ! Use internal read for robust string-to-integer conversion
      read(line_data(year_idx), *, iostat=ierr) year_int
      if (ierr /= 0) call error_stop("Error reading 'year' as integer")
      
      read(line_data(month_idx), *, iostat=ierr) month_int
      if (ierr /= 0) call error_stop("Error reading 'month' as integer")
      
      read(line_data(day_idx), *, iostat=ierr) day_int
      if (ierr /= 0) call error_stop("Error reading 'day' as integer")

      read(line_data(hour_idx), *, iostat=ierr) hour_int
      if (ierr /= 0) call error_stop("Error reading 'hour' as integer")

      read(line_data(minute_idx), *, iostat=ierr) minute_int
      if (ierr /= 0) call error_stop("Error reading 'minute' as integer")

      ! 5. Format integers to strings as %02d (and %04d for year)
      ! Use internal write for robust integer-to-string formatting
      write(year_str, "(I4.4)") year_int
      write(month_str, "(I2.2)") month_int
      write(day_str, "(I2.2)") day_int
      write(hour_str, "(I2.2)") hour_int
      write(minute_str, "(I2.2)") minute_int
      
      ! 6. Construct the final filename
      output_filename = trim(year_str) // trim(month_str) // trim(day_str) &
                      & // "_" // trim(hour_str) // trim(minute_str) // ".nc"
                      
      ! --- Process the data (example: print it) ---
      print '("--- Processing Row ---")'
      print '("  Time: ", A, ", Gas: ", A, ", Path: ", A)', trim(time_str), trim(gas_str), trim(path_str)
      print '("  Filename: ", A)', trim(output_filename)
      !
      ! You would add logic here to use 'output_filename' and 'path_str' to read the netCDF file
      !
      
    end do

    ! Close the file
    call close_csv(csv)
    
  end subroutine example_csv_reader
  
end module io_manager