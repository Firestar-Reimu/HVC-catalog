```text
Rotation curve C extension module
=================================

Files:
- rotation_curve.c   : C source for Python extension module "rotation_curve_c"
- setup.py           : build script using setuptools
- example_usage.py   : example of how to import and call calc_v_dev from Python

Build:
    python setup.py build_ext --inplace

This will generate a shared library like rotation_curve_c*.so (Linux/macOS) or rotation_curve_c*.pyd (Windows)
in the current directory. Then you can import it:

Usage example (see example_usage.py):
    import rotation_curve_c
    vmax, vmin = rotation_curve_c.calc_v_dev(230, 10, model='univ')
    print(vmax, vmin)

The function signature available from Python:
    calc_v_dev(l, b, h=5.0, r_gal=20.0, r_sun=8.5, n_sample=1000, model='univ', dev=50.0)
Return:
    (v_max, v_min) as floats

Notes:
- Angles l and b are in degrees.
- model choices: 'simple', 'univ', 'linear', 'linear2'.
- Implementation mirrors the original Python code's numeric behavior (within C floating point limits).
```