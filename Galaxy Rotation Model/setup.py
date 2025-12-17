from setuptools import setup, Extension
import numpy as np

module = Extension(
    'rotation_model_c',
    sources=['rotation_model_c.c'],
    extra_compile_args=['-O3'],
)

setup(
    name='rotation_model_c',
    description='C extension of calc_v_dev from rotation_model_py.py',
    ext_modules=[module],
)

module = Extension(
    'process_cube_c',
    sources=['process_cube_c.c'],
    extra_compile_args=['-O3'],
    include_dirs=[np.get_include()],
)

setup(
    name='process_cube_c',
    description='C extension for fast FITS cube masking based on deviation velocity',
    ext_modules=[module],
)