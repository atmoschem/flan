program test_write
  use io_manager, only: write_1d_text
  use linear_algebra, only: diag
  use iso_fortran_env, only: wp => real64
  implicit none
  
  real(wp), dimension(5) :: test_data, r_diag
  real(wp), dimension(5,5) :: r_mat
  character(len=256) :: test_file = "test_hsp.txt"
  real(wp) :: obs_err_sd = 1.0_wp
  real(wp) :: model_err_sd = 2.0_wp
  
  test_data = [1.1_wp, 2.2_wp, 3.3_wp, 4.4_wp, 5.5_wp]
  
  print *, "Testing write_1d_text..."
  call write_1d_text(test_file, test_data)
  
  print *, "Testing R matrix calculation..."
  r_diag = obs_err_sd**2 + model_err_sd**2
  r_mat = diag(r_diag)
  
  print *, "R diagonal: ", r_diag(1)
  print *, "R shape: ", shape(r_mat)
  print *, "R(1,1): ", r_mat(1,1)
  print *, "R(1,2): ", r_mat(1,2)
  
  if (abs(r_mat(1,1) - 5.0_wp) < 1e-10_wp .and. abs(r_mat(1,2)) < 1e-10_wp) then
    print *, "SUCCESS: R matrix correctly created."
  else
    print *, "FAILURE: R matrix error."
    stop 1
  end if
  
end program test_write
