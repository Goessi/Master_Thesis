from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

setup(
    ext_modules = cythonize(Extension(
    'Quaternion',
    sources=['Quaternion.pyx'],
    language='c++',
    include_dirs=[numpy.get_include()]
))
)
