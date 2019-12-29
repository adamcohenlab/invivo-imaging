from setuptools import setup, Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

import numpy

setup(
    name="l1_tf_C",
    ext_modules=cythonize([Extension("l1_tf_C.c_l1_tf", ["l1_tf_C/c_l1_tf.pyx"],
                                     include_dirs=[numpy.get_include()],
                                     libraries=["blas", "lapack", "m"]),]),
    cmdclass={"build_ext": build_ext},
    packages=["l1_tf_C",],
)
