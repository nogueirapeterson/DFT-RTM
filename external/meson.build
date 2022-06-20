pyplot_fortran_lib = library('pyplot_fortran',
  'pyplot-fortran/src/pyplot_module.F90',
  fortran_args : [ '-cpp', '-DREAL32' ])

pyplot_fortran_dep = declare_dependency(link_with: pyplot_fortran_lib)
