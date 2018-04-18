! some wrapper functions used for interoperability with C/C++

FUNCTION get_string(c_rc) bind(c, name='get_string')
  ! ref: https://stackoverflow.com/a/44424572
  USE, INTRINSIC :: iso_c_binding
  IMPLICIT NONE
  INTEGER(c_int), INTENT(out) :: c_rc           ! <- Pass by reference; acts as return code in this example.
  TYPE(c_ptr) :: get_string                     ! <- C_PTR to pass back to C
  CHARACTER(len=:), ALLOCATABLE, TARGET, SAVE :: fortstring   ! <- Allocatable/any length is fine

  fortstring = "Append C_NULL_CHAR to any Fortran string constant."//C_NULL_CHAR
  get_string = c_loc(fortstring)                ! <- C_LOC intrinsic gets loc of our string.
  c_rc = 1                                      ! <- Set the return code value.
END FUNCTION get_string


INTEGER(kind=c_int) FUNCTION put_string(instr) bind(C,name='put_string')
  ! ref: https://stackoverflow.com/a/5609853
  USE, INTRINSIC :: iso_c_binding
  CHARACTER(kind=c_char), DIMENSION(*), INTENT(IN) :: instr
  INTEGER :: len
  INTEGER :: i

  len=0
  DO
     IF (instr(len+1) == C_NULL_CHAR) EXIT
     len = len + 1
  END DO


  PRINT *, 'In Fortran:'
  PRINT *, 'Got string: ', (instr(i),i=1,len)
  put_string = len
END FUNCTION put_string


FUNCTION get_name_c(pos,c_rc) bind(c, name='get_name_c')
  USE, INTRINSIC :: iso_c_binding
  USE SUEWS_Driver,ONLY:get_name
  IMPLICIT NONE
  INTEGER(c_int), INTENT(in) :: pos
  INTEGER(c_int), INTENT(out) :: c_rc           ! <- Pass by reference; acts as return code in this example.
  TYPE(c_ptr) :: get_name_c                     ! <- C_PTR to pass back to C
  CHARACTER(len=:), ALLOCATABLE, TARGET, SAVE :: fortstring   ! <- Allocatable/any length is fine
  CHARACTER(len=15)::name


  CALL get_name(pos,name) ! <- get variable name from varList
  fortstring = TRIM(name)//C_NULL_CHAR
  get_name_c = c_loc(fortstring)                ! <- C_LOC intrinsic gets loc of our string.
  c_rc = 1                                      ! <- Set the return code value.
END FUNCTION get_name_c
