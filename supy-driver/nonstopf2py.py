"""
TS 19 Feb 2019 :
adopted this script below to resolve the Fortran `stop` issue:
Credit to: https://github.com/joezuntz/pycamb/blob/master/nonstopf2py.py
This is an adapted version with added modification for windows,
where several GNU C libs are missing that require windows alternatives.

original header:

   An non-stop f2py. When the fortran program calls stop, it is trapped with a long jmp,
   and the error is converted to a python error.

   A lot of inspiration came from patch_f2py.py in James Kermode's quippy
     (http://www.jrkermode.co.uk/quippy)

   currently supported fortran runtimes are:
    ifcore (intel)
    gfortran (gnu)

   For new runtimes, compile a small f90 that uses stop

   boo.f90:
      program boo
          stop "boo the dog stops here"
      end program boo

   and

     objdump -t a.o

   to find the function to override.

  -- Yu Feng <yfeng1@cmu.edu> @ McWilliam Center, Carnegie Mellon 2012

"""

import platform

from numpy import f2py
print('Customised f2py loaded!')


# trap the fortran STOP methods.
# luckily they are not builtin/inlined.
# we do not need to free 'message'. it appears to be staticly allocated
# by the compiler.

f2py.rules.module_rules['modulebody'] = f2py.rules.module_rules['modulebody'].replace(
    '#includes0#\n',
    r"""#includes0#
#include <setjmp.h>
static char * _error;
static jmp_buf _env;
void for_stop_core(char * message, int len) {
  _error = strndup(message, len);
  longjmp(_env, 1);
}
void _gfortran_stop_string(char * message, int len) {
  _error = strndup(message, len);
  longjmp(_env, 1);
}
     """)

# here we fight the leak as f2py will no longer always return.
# the easiest way is to first construct the return tuple,
# then free them all
f2py.rules.routine_rules['body'] = f2py.rules.routine_rules['body'].replace(
    """\t\tif (f2py_success) {
#pyobjfrom#
/*end of pyobjfrom*/
\t\tCFUNCSMESS(\"Building return value.\\n\");
\t\tcapi_buildvalue = Py_BuildValue(\"#returnformat#\"#return#);
/*closepyobjfrom*/
#closepyobjfrom#
\t\t} /*if (f2py_success) after callfortranroutine*/""",

    """\t\t{
#pyobjfrom#
/*end of pyobjfrom*/
\t\tCFUNCSMESS(\"Building return value.\\n\");
\t\tcapi_buildvalue = Py_BuildValue(\"#returnformat#\"#return#);
/*closepyobjfrom*/
#closepyobjfrom#
\t\tif(!f2py_success) {
\t\t\tPy_XDECREF(capi_buildvalue);
\t\t\tcapi_buildvalue = NULL;
\t\t}
\t\t}
/*if (f2py_success) after callfortranroutine*/
"""
)

# the actual function call. free _error as PyErr_SetString will copy it.
f2py.rules.routine_rules['body'] = f2py.rules.routine_rules['body'].replace(
    '#callfortranroutine#\n',
    r"""
       if(setjmp(_env)) {
         PyErr_SetString(PyExc_RuntimeError, _error);
         free(_error);
       } else {
         #callfortranroutine#
       }
   """)

# distinguish platform to handle lib missing issues on windows
sysname = platform.system()

# change `setjmp` and `longjmp` to windows compliant versions
# http://www.agardner.me/golang/windows/cgo/64-bit/setjmp/longjmp/2016/02/29/go-windows-setjmp-x86.html
if sysname == 'Windows':
    f2py.rules.routine_rules['body'] = f2py.rules.routine_rules['body'].replace(
        r'setjmp(_env)',
        r'__builtin_setjmp(_env)'
    )
    f2py.rules.module_rules['modulebody'] = f2py.rules.module_rules['modulebody'].replace(
        r'longjmp(_env)',
        r'__builtin_longjmp(_env)'
    )

    # add an implementation of `strndup`
    # https://github.com/noahp/cflow-mingw/blob/4e3f48c6636f4e54e2f30671eaeb3c8bc96fd4a4/cflow-1.4/gnu/strndup.c
    f2py.rules.module_rules['modulebody'] = f2py.rules.module_rules['modulebody'].replace(
        r'#include <setjmp.h>',
        r"""
#include <setjmp.h>
#include <stdlib.h>
char *
strndup (char const *s, size_t n)
{
  size_t len = strnlen (s, n);
  char *new = malloc (len + 1);

  if (new == NULL)
    return NULL;

  new[len] = '\0';
  return memcpy (new, s, len);
}
"""
    )

