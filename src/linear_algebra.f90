
module linear_algebra
  use iso_fortran_env, only: wp => real64
 ! CORRECT: Only import 'inv'. 'norm2' is a built-in Fortran function.
  use stdlib_linalg, only: inv, diag
  implicit none
  private
  public :: wp, vec_norm, invert_matrix, diag

contains

  function vec_norm(x) result(norm)
    real(wp), intent(in) :: x(:)
    real(wp) :: norm
   ! This will now correctly use the intrinsic norm2 function from the Fortran language itself.
    norm = norm2(x)
  end function vec_norm

  subroutine invert_matrix(A)
    real(wp), intent(inout) :: A(:,:)
   ! This correctly uses the high-level matrix inversion function from stdlib.
    A = inv(A)
  end subroutine invert_matrix

end module linear_algebra

