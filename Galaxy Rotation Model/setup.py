from setuptools import setup, Extension

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