module flim
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, flim!"
  end subroutine say_hello
end module flim
