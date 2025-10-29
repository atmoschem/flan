module io_manager
  use iso_fortran_env, only: wp => real64
  use netcdf
  use stdlib_string_type, only: string_type
  use stdlib_strings, only: to_string !, strip
  use stdlib_ascii, only: to_lower
  use fpm_strings, only: split    ! 'split' is in fpm_strings
  use stdlib_error, only: error_stop
  implicit none
  private
 ! CORRECT: Make both subroutines public
  public :: write_matrix_netcdf, read_matrix_netcdf

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

subroutine read_csv_obs(filename, path_col, time_col, gas_col, bg_col, &
                          year_col, month_col, day_col, hour_col, minute_col)
    ! Inputs
    character(len=*), intent(in) :: filename
    ! Outputs
    character(len=:), allocatable, intent(out) :: path_col(:)
    character(len=:), allocatable, intent(out) :: time_col(:)
    real(wp), allocatable, intent(out) :: gas_col(:)
    real(wp), allocatable, intent(out) :: bg_col(:)
    ! Date/time columns are now formatted characters
    character(len=:), allocatable, intent(out) :: year_col(:)
    character(len=:), allocatable, intent(out) :: month_col(:)
    character(len=:), allocatable, intent(out) :: day_col(:)
    character(len=:), allocatable, intent(out) :: hour_col(:)
    character(len=:), allocatable, intent(out) :: minute_col(:)

    ! Local variables
    integer :: unit_num, io_stat
    character(len=2048) :: line_buffer
    type(string_type) :: header_line, data_line
    character(len=:), allocatable :: headers(:)
    character(len=:), allocatable :: data_values(:)
    integer :: i, num_cols
    integer :: path_idx, time_idx, gas_idx, bg_idx, year_idx, &
               month_idx, day_idx, hour_idx, minute_idx
    
    ! Temporary variables for reading and formatting
    integer :: int_year, int_month, int_day, int_hour, int_minute
    character(len=4) :: char_year
    character(len=2) :: char_month, char_day, char_hour, char_minute

    ! Temporary allocatable arrays for building output
    character(len=512), allocatable :: path_tmp(:), time_tmp(:)
    real(wp), allocatable :: gas_tmp(:), bg_tmp(:)
    character(len=4), allocatable :: year_tmp(:)
    character(len=2), allocatable :: month_tmp(:), day_tmp(:), &
                                     hour_tmp(:), minute_tmp(:)

    ! Initialize indices to -1 (not found)
    path_idx = -1; time_idx = -1; gas_idx = -1; bg_idx = -1; year_idx = -1; &
    month_idx = -1; day_idx = -1; hour_idx = -1; minute_idx = -1

    ! Allocate temporary arrays
    allocate(path_tmp(0)); allocate(time_tmp(0))
    allocate(gas_tmp(0)); allocate(bg_tmp(0))
    allocate(year_tmp(0)); allocate(month_tmp(0)); allocate(day_tmp(0)); &
    allocate(hour_tmp(0)); allocate(minute_tmp(0))

    open(newunit=unit_num, file=filename, action="read", status="old", iostat=io_stat)
    if (io_stat /= 0) call error_stop("Error: Could not open file: " // trim(filename))

    ! 1. Read header row
    read(unit_num, '(a)', iostat=io_stat) line_buffer
    if (io_stat /= 0) call error_stop("Error: Could not read header from file: " // trim(filename))
    
    header_line = to_string(strip(line_buffer))
    call split(',', headers)
    num_cols = size(headers)

    ! 2. Find column indices (case-insensitive)
    do i = 1, num_cols
      select case (to_lower(strip(headers(i))))
        case ('path_col')
          path_idx = i
        case ('time_col')
          time_idx = i
        case ('gas_col')
          gas_idx = i
        case ('bg_col')
          bg_idx = i
        case ('year')
          year_idx = i
        case ('month')
          month_idx = i
        case ('day')
          day_idx = i
        case ('hour')
          hour_idx = i
        case ('minute')
          minute_idx = i
      end select
    end do

    ! 3. Check if all required columns were found
    if (path_idx == -1) call error_stop("Error: 'path_col' not found in " // trim(filename))
    if (time_idx == -1) call error_stop("Error: 'time_col' not found in " // trim(filename))
    if (gas_idx == -1) call error_stop("Error: 'gas_col' not found in " // trim(filename))
    if (bg_idx == -1) call error_stop("Error: 'bg_col' not found in " // trim(filename))
    if (year_idx == -1) call error_stop("Error: 'year' not found in " // trim(filename))
    if (month_idx == -1) call error_stop("Error: 'month' not found in " // trim(filename))
    if (day_idx == -1) call error_stop("Error: 'day' not found in " // trim(filename))
    if (hour_idx == -1) call error_stop("Error: 'hour' not found in " // trim(filename))
    if (minute_idx == -1) call error_stop("Error: 'minute' not found in " // trim(filename))

    ! 4. Read data rows
    data_loop: do
      read(unit_num, '(a)', iostat=io_stat) line_buffer
      if (io_stat /= 0) exit data_loop ! Exit on end-of-file or read error

      data_line = strip(line_buffer)
      if (len(data_line) == 0) cycle data_loop ! Skip empty lines
      
      call split(',', data_values)

      if (size(data_values) == num_cols) then
        ! Append string data
        path_tmp = [path_tmp, strip(data_values(path_idx))]
        time_tmp = [time_tmp, strip(data_values(time_idx))]
        
        ! Append real data
        gas_tmp = [gas_tmp, 0.0_wp] ! Allocate space
        read(data_values(gas_idx), *) gas_tmp(size(gas_tmp))
        
        bg_tmp = [bg_tmp, 0.0_wp] ! Allocate space
        read(data_values(bg_idx), *) bg_tmp(size(bg_tmp))

        ! Read integer date/time data
        read(data_values(year_idx), *) int_year
        read(data_values(month_idx), *) int_month
        read(data_values(day_idx), *) int_day
        read(data_values(hour_idx), *) int_hour
        read(data_values(minute_idx), *) int_minute

        ! Format as zero-padded strings and append
        write(char_year, '(i4.4)') int_year
        year_tmp = [year_tmp, char_year]
        
        write(char_month, '(i2.2)') int_month
        month_tmp = [month_tmp, char_month]

        write(char_day, '(i2.2)') int_day
        day_tmp = [day_tmp, char_day]

        write(char_hour, '(i2.2)') int_hour
        hour_tmp = [hour_tmp, char_hour]

        write(char_minute, '(i2.2)') int_minute
        minute_tmp = [minute_tmp, char_minute]
        
      else
         print *, "Warning: Skipping malformed line: ", trim(line_buffer)
      end if
    end do data_loop

    close(unit_num)

    ! 5. Move data from temporary arrays to final output arrays
    path_col = path_tmp
    time_col = time_tmp
    gas_col = gas_tmp
    bg_col = bg_tmp
    year_col = year_tmp
    month_col = month_tmp
    day_col = day_tmp
    hour_col = hour_tmp
    minute_col = minute_tmp

  end subroutine read_csv_obs

  
end module io_manager