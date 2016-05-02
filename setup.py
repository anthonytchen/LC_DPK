# setup.py

from distutils.core import setup, Extension

setup(name="Chemical",
      py_modules=['Chemical'], 
      ext_modules=[Extension("_Chemical",
                             ["Chemical.i","Chemical.cpp"],
                             swig_opts=['-c++'],
                         )
               ]
)

#setup(name="Skin",
#      py_modules=['Skin'], 
#      ext_modules=[Extension("_Skin",
#                             ["Skin.i", "Chemical.cpp", "Grid.cpp", "StraCorn.cpp", "ViaEpd.cpp", "Dermis.cpp", "Blood.cpp", "Skin.cpp"],
#                             swig_opts=['-c++'],
#                             library_dirs=['/usr/lib', '/usr/local/lib'],
#                            libraries=['gsl', 'gslcblas', 'sundials_cvode', 'sundials_nvecserial', 'm'],
#                             include_dirs=['/usr/include', '/usr/local/include'],
#                         )
#               ]
#)

setup(name="Skin_S3VDB",
      py_modules=['Skin_S3VDB'], 
      ext_modules=[Extension("_Skin_S3VDB",
                             ["Skin_S3VDB.i","Chemical.cpp", "except.cpp", "Grid.cpp", "Comp.cpp", "Vehicle.cpp", "Sebum.cpp", "SurSebum.cpp", "StraCorn.cpp", "ViaEpd.cpp", "Dermis.cpp", "Skin.cpp", "Skin_S3VDB.cpp", "Blood.cpp"],
                             swig_opts=['-c++'],
                             library_dirs=['/usr/lib', '/usr/local/lib'],
                             libraries=['sundials_cvode', 'sundials_nvecserial', 'm'],
                             include_dirs=['/usr/include', '/usr/local/include'],
                         )
               ]
)
