module io_manager
  use iso_fortran_env, only: wp => real64
  use netcdf
  use csv_module, only: csv_file, open_csv, close_csv, read_line
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

! --- Example CSV Reader using csv-fortran ---
  subroutine example_csv_reader(filename)
    character(len=*), intent(in) :: filename
    type(csv_file) :: csv
    character(len=:), allocatable :: line_data(:)
    logical :: at_end
    integer :: i

    ! Open the CSV file
    call csv%open(csv, filename, 'read')
    if (.not. csv%is_open()) then
      call error_stop("Error: Could not open CSV file: " // trim(filename))
    end if

    ! Read header (optional, but good practice)
    call read_line(csv, line_data, at_end)
    if (at_end) then
      print *, "Warning: CSV file is empty."
      call close_csv(csv)
      return
    end if
    print *, "CSV Headers:"
    print '("  ", A)', line_data

    ! Loop through data rows
    do
      call read_line(csv, line_data, at_end)
      if (at_end) exit ! Exit loop at end of file
      
      ! Process the line_data array (which contains strings for each column)
      ! Example: Print the first 3 columns
      print *, "Data Row: "
      do i = 1, min(size(line_data), 3)
         print '("    Col", I0, ": ", A)', i, trim(line_data(i))
      end do
      
      ! You would add logic here to read(line_data(i), *) into real/integer variables
      
    end do

    ! Close the file
    call close_csv(csv)
    
  end subroutine example_csv_reader
  
end module io_manager