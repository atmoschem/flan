module io_manager
  use iso_fortran_env, only: wp => real64
  use netcdf
  use csv_module
  use stdlib_error, only: error_stop
  use datetime_module, only:datetime
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
  subroutine example_csv_reader(csvname)
    character(len=*), intent(in) :: csvname
    type(csv_file) :: csv
    character(len=30),dimension(:),allocatable :: header
    logical :: at_end
    integer :: i, ierr
    
    ! type time
    type(datetime) :: a

    !Variables to hold data from a row
    character(len=100), dimension(:), allocatable :: path
    real(wp),dimension(:),allocatable :: gas, bg
    integer,dimension(:),allocatable :: year, month, day, hour, minute
    character(len=100) :: output_filename
    integer, dimension(:), allocatable :: itypes


    call csv%read(csvname,header_row=1,status_ok=at_end)

    ! get the header and type info
    call csv%get_header(header,at_end)
    call csv%variable_types(itypes,at_end)

    call csv%get(1,year,at_end)
    call csv%get(2,month,at_end)
    call csv%get(3,day,at_end)
    call csv%get(4,hour,at_end)
    call csv%get(5,minute,at_end)
    call csv%get(6,path,at_end)
    call csv%get(7,gas,at_end)
    call csv%get(8,bg,at_end)

    call csv%destroy()
     
    a = datetime(year(1), month(1), day(1), hour(1), minute(1), 0) 

    print *, a%isoformat() 
    ! 6. Construct the final filename
!      output_filename = trim(year_str) // trim(month_str) // trim(day_str) &
!                      & // "_" // trim(hour_str) // trim(minute_str) // ".nc"
                      
    end subroutine example_csv_reader

end module io_manager