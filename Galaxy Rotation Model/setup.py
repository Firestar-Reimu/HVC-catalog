from setuptools import setup, Extension

module = Extension(
    'rotation_curve_c',
    sources=['rotation_curve.c'],
    extra_compile_args=[],
)

setup(
    name='rotation_curve_c',
    version='0.1',
    description='C extension of calc_v_dev from rotation_curve.py',
    ext_modules=[module],
)

module = Extension(
    'process_cube_c',
    sources=['process_cube_c.c'],
    extra_compile_args=['-O3'],
)

setup(
    name='process_cube_c',
    version='0.1',
    description='C extension for fast FITS cube masking based on deviation velocity',
    ext_modules=[module],
)