# setup.py

from distutils.core import setup, Extension

setup(name="Chemical",
      py_modules=['Chemical'], 
      ext_modules=[Extension("_Chemical",
                             ["Chemical.i","Chemical.cpp"],
                             swig_opts=['-c++'],
							 library_dirs=['/usr/lib', '/usr/local/lib', '/user/HS104/tc0008/work/lib/'],
                             libraries=['sundials_cvode', 'sundials_nvecserial', 'm'],
							 include_dirs=['/usr/include', '/usr/local/include', '/user/HS104/tc0008/work/include/'],
                         )
               ]
)

setup(name="Config",
      py_modules=['Config'], 
      ext_modules=[Extension("_Config",
                             ["Config.i","Config.cpp"],
                             swig_opts=['-c++'],
							 library_dirs=['/usr/lib', '/usr/local/lib', '/user/HS104/tc0008/work/lib/'],
                             libraries=['sundials_cvode', 'sundials_nvecserial', 'm'],
							 include_dirs=['/usr/include', '/usr/local/include', '/user/HS104/tc0008/work/include/'],
                         )
               ]
)


setup(name="Skin_Setup",
      py_modules=['Skin_Setup'], 
      ext_modules=[Extension("_Skin_Setup",
                             ["Skin_Setup.i",
							 "except.cpp", "Config.cpp", "Chemical.cpp", "Grid.cpp", "Comp.cpp",  "Vehicle.cpp",  "Sebum.cpp",  "SurSebum.cpp",  "StraCorn.cpp",  "ViaEpd.cpp",  "Dermis.cpp",  "Skin.cpp",  "Skin_Setup.cpp",  "Blood.cpp"],
                             swig_opts=['-c++'],
                             library_dirs=['/usr/lib', '/usr/local/lib', '/user/HS104/tc0008/work/lib/','/usr/local/lib/python2.7/dist-packages/numpy/core/lib'],
                             libraries=['sundials_cvode', 'sundials_nvecserial', 'm'],
                             include_dirs=['/usr/include', '/usr/local/include', '/user/HS104/tc0008/work/include/', '/usr/local/lib/python2.7/dist-packages/numpy/core/include'],
                         )
               ]
)
