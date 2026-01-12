module io_manager
  use iso_fortran_env, only: wp => real64
  use netcdf
  use csv_module
  use stdlib_error, only: error_stop
  use datetime_module, only:datetime
  implicit none

  ! CORRECT: Make both subroutines public
  public :: write_1d_netcdf, write_2d_netcdf, write_3d_netcdf, read_1d_netcdf, read_2d_netcdf, read_3d_netcdf
  public :: example_csv_reader
  public :: logical_to_string, write_1d_text

contains

  subroutine check(status)
   ! Helper subroutine for error checking
    integer, intent(in) :: status
    if (status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "NetCDF Error"
    end if
  end subroutine check


  subroutine write_1d_netcdf(filename, data_matrix, varname)
    character(len=*), intent(in) :: filename
    real(wp),  intent(in) :: data_matrix(:)
    character(len=*), intent(in) :: varname

    integer :: ncid, dim1_id, var_id
    integer :: status
    integer :: dim1_size

    dim1_size = size(data_matrix, 1)

    status = nf90_create(filename, nf90_clobber, ncid)
    call check(status)

    status = nf90_def_dim(ncid, "dim1", dim1_size, dim1_id)
    call check(status)

    status = nf90_def_var(ncid, varname, nf90_double, [dim1_id], var_id)
    call check(status)

    status = nf90_enddef(ncid)
    call check(status)

    status = nf90_put_var(ncid, var_id, data_matrix)
    call check(status)

    status = nf90_close(ncid)
    call check(status)
  end subroutine write_1d_netcdf


  subroutine write_2d_netcdf(filename, data_matrix, varname)
    character(len=*), intent(in) :: filename
    real(wp),  intent(in) :: data_matrix(:,:)
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

    status = nf90_def_var(ncid, varname, nf90_double, [dim1_id, dim2_id], var_id)
    call check(status)

    status = nf90_enddef(ncid)
    call check(status)

    status = nf90_put_var(ncid, var_id, data_matrix)
    call check(status)

    status = nf90_close(ncid)
    call check(status)
  end subroutine write_2d_netcdf


  subroutine write_3d_netcdf(filename, data_matrix, varname)
    character(len=*), intent(in) :: filename
    real(wp),  intent(in) :: data_matrix(:,:,:)
    character(len=*), intent(in) :: varname

    integer :: ncid, dim1_id, dim2_id, dim3_id, var_id
    integer :: status
    integer :: dim1_size, dim2_size, dim3_size

    dim1_size = size(data_matrix, 1)
    dim2_size = size(data_matrix, 2)
    dim3_size = size(data_matrix, 3)

    status = nf90_create(filename, nf90_clobber, ncid)
    call check(status)

    status = nf90_def_dim(ncid, "dim1", dim1_size, dim1_id)
    call check(status)
    status = nf90_def_dim(ncid, "dim2", dim2_size, dim2_id)
    call check(status)
    status = nf90_def_dim(ncid, "dim3", dim3_size, dim3_id)
    call check(status)

    status = nf90_def_var(ncid, varname, nf90_double, [dim1_id, dim2_id, dim3_id], var_id)
    call check(status)

    status = nf90_enddef(ncid)
    call check(status)

    status = nf90_put_var(ncid, var_id, data_matrix)
    call check(status)

    status = nf90_close(ncid)
    call check(status)
  end subroutine write_3d_netcdf


  subroutine read_1d_netcdf(filename, data_matrix, varname)
    character(len=*), intent(in) :: filename
    real(wp), allocatable, intent(inout) :: data_matrix(:)
    character(len=*), intent(in) :: varname

    integer :: ncid, var_id, status
    integer :: dimids(1), dim_size
    integer :: ndims ! Added for dimension check

   ! 1. Open the file for reading.
    status = nf90_open(filename, nf90_nowrite, ncid)
    call check(status)

   ! 2. Get the ID of the variable by its name.
    status = nf90_inq_varid(ncid, varname, var_id)
    call check(status)

    ! Check number of dimensions
    status = nf90_inquire_variable(ncid, var_id, ndims=ndims)
    call check(status)
    
    if (ndims /= 1) then
        print *, "Error: Variable ", trim(varname), " has ", ndims, " dimensions, expected 1"
        stop "Dimension mismatch"
    end if

    ! Get dimensions from file
    status = nf90_inquire_variable(ncid, var_id, dimids=dimids)
    call check(status)
    status = nf90_inquire_dimension(ncid, dimids(1), len=dim_size)
    call check(status)

    ! Check if reallocation is needed
    if (allocated(data_matrix)) then
       if (size(data_matrix) /= dim_size) then
           deallocate(data_matrix)
       end if
    end if

    ! Allocate if not allocated
    if (.not. allocated(data_matrix)) then
      allocate(data_matrix(dim_size))
    end if

   ! 3. Read the data from the variable into the array.
    status = nf90_get_var(ncid, var_id, data_matrix)
    call check(status)

   ! 4. Close the file.
    status = nf90_close(ncid)
    call check(status)
  end subroutine read_1d_netcdf


  subroutine read_2d_netcdf(filename, data_matrix, varname)
    character(len=*), intent(in) :: filename
    real(wp),  intent(inout) :: data_matrix(:,:)
    character(len=*), intent(in) :: varname

    integer :: ncid, var_id, status

    status = nf90_open(filename, nf90_nowrite, ncid)
    call check(status)

    status = nf90_inq_varid(ncid, varname, var_id)
    call check(status)

    status = nf90_get_var(ncid, var_id, data_matrix)
    call check(status)

    status = nf90_close(ncid)
    call check(status)
  end subroutine read_2d_netcdf


  subroutine read_3d_netcdf(filename, data_matrix, varname)
    character(len=*), intent(in) :: filename
    real(wp), allocatable, intent(inout) :: data_matrix(:,:,:)
    character(len=*), intent(in) :: varname

    integer :: ncid, var_id, status
    integer :: dimids(3), dim_sizes(3)
    integer :: ndims

   ! 1. Open the file for reading.
    status = nf90_open(filename, nf90_nowrite, ncid)
    call check(status)

   ! 2. Get the ID of the variable by its name.
    status = nf90_inq_varid(ncid, varname, var_id)
    call check(status)

    ! Check number of dimensions
    status = nf90_inquire_variable(ncid, var_id, ndims=ndims)
    call check(status)
    
    if (ndims /= 3) then
        print *, "Error: Variable ", trim(varname), " has ", ndims, " dimensions, expected 3"
        stop "Dimension mismatch"
    end if

    ! Get dimensions from file
    status = nf90_inquire_variable(ncid, var_id, dimids=dimids)
    call check(status)
    status = nf90_inquire_dimension(ncid, dimids(1), len=dim_sizes(1))
    call check(status)
    status = nf90_inquire_dimension(ncid, dimids(2), len=dim_sizes(2))
    call check(status)
    status = nf90_inquire_dimension(ncid, dimids(3), len=dim_sizes(3))
    call check(status)

    ! Check if reallocation is needed
    if (allocated(data_matrix)) then
       if (size(data_matrix, 1) /= dim_sizes(1) .or. &
           size(data_matrix, 2) /= dim_sizes(2) .or. &
           size(data_matrix, 3) /= dim_sizes(3)) then
           deallocate(data_matrix)
       end if
    end if

    ! Allocate if not allocated
    if (.not. allocated(data_matrix)) then
      allocate(data_matrix(dim_sizes(1), dim_sizes(2), dim_sizes(3)))
    end if

   ! 3. Read the data from the variable into the array.
    status = nf90_get_var(ncid, var_id, data_matrix)
    call check(status)

    status = nf90_close(ncid)
    call check(status)
  end subroutine read_3d_netcdf


! --- Example CSV Reader using csv-fortran ---
subroutine example_csv_reader(csvname, year, month, day, hour, minute, second, path, gas, bg, err)
    character(len=*), intent(in) :: csvname

! These are now output arguments
    integer, dimension(:), allocatable, intent(out) :: year, month, day, hour, minute, second
    character(len=256), dimension(:), allocatable, intent(out) :: path
    real(wp), dimension(:), allocatable, intent(out) :: gas, bg, err

    type(csv_file) :: csv
    character(len=30),dimension(:),allocatable :: header
    logical :: at_end
    integer :: i, ierr, num_records

    ! time
    type(datetime) :: a

    character(len=100) :: output_filename
    integer, dimension(:), allocatable :: itypes

    print *, "Reading CSV file: ", trim(csvname)

    call csv%read(csvname,header_row=1,status_ok=at_end)

    ! get the header and type info
    call csv%get_header(header,at_end)
    call csv%variable_types(itypes,at_end)

    call csv%get(1,year,at_end)
    call csv%get(2,month,at_end)
    call csv%get(3,day,at_end)
    call csv%get(4,hour,at_end)
    call csv%get(5,minute,at_end)
    ! second is not in CSV
    call csv%get(6,path,at_end)
    call csv%get(7,gas,at_end)
    call csv%get(8,bg,at_end)
    call csv%get(9,err,at_end) 

    ! Get the number of records read
    num_records = size(year)

    ! Allocate and initialize second
    allocate(second(num_records))
    second = 0


    call csv%destroy()

    ! This part just prints the datetime of the first record
    if (num_records > 0) then
        print *, "--- CSV Summary ---"
        print *, "Total records: ", num_records
        print *, "First Record:"
        print *, "  Time: ", year(1), month(1), day(1), hour(1), minute(1), second(1)
        print *, "  Gas: ", gas(1)
        print *, "  Bg: ", bg(1)
        print *, "  Path: ", path(1)

        print *, "Last Record:"
        print *, "  Time: ", year(num_records), month(num_records), day(num_records), &
                             hour(num_records), minute(num_records), second(num_records)
        print *, "  Gas: ", gas(num_records)
        print *, "  Bg: ", bg(num_records)
        print *, "  Path: ", path(num_records)

        a = datetime(year(1), month(1), day(1), hour(1), minute(1), 0)
        print *, "First receptor time object: ", a%isoformat()
    end if
     
    ! 6. Construct the final filename
!    output_filename = trim(year_str) // trim(month_str) // trim(day_str) &
!                      // "_" // trim(hour_str) // trim(minute_str) // ".nc"

    end subroutine example_csv_reader
    
    function logical_to_string(l) result(s)
      logical, intent(in) :: l
      character(len=5) :: s
      if (l) then
        s = "true"
      else
        s = "false"
      end if
    end function logical_to_string

    subroutine write_1d_text(filename, data_matrix)
      character(len=*), intent(in) :: filename
      real(wp), intent(in) :: data_matrix(:)
      integer :: unit_num, i, ios

      open(newunit=unit_num, file=filename, status='replace', action='write', iostat=ios)
      if (ios /= 0) then
        print *, "Error opening file for writing: ", trim(filename)
        return
      end if

      do i = 1, size(data_matrix)
        write(unit_num, '(F15.5)') data_matrix(i)
      end do

      close(unit_num)
    end subroutine write_1d_text

end module io_manager