project('DFT-RTM', 'fortran',
  version : '0.1',
  default_options : ['warning_level=3'])

subdir('external')

dft_lib = library('dft',
           'src/mdle_adds.f90',
           'src/mdle_io_utils.f90')

dft_dep = declare_dependency(link_with : dft_lib)

executable('DFT-RTM',
           'app/main.f90',
           dependencies : dft_dep,
           install : true)

test1 = executable('DFT-RTM-test',
           'test/test.f90',
           dependencies : [ pyplot_fortran_dep, dft_dep ])

test('testing 51x51 model', test1)
