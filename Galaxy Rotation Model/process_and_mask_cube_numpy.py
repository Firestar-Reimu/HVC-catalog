#!/usr/bin/env python3
"""
process_and_mask_cube.py

Simplified processing script:

- Uses spectral_cube exclusively to extract the spectral axis in km/s:
    from spectral_cube import SpectralCube as sc
    cube = sc.read(file).with_spectral_unit(u.km / u.s)
    velocities = cube.spectral_axis.value

- For each spatial pixel (x,y) compute Galactic coordinates (l,b) using the cube WCS,
  call calc_v_dev(l,b,...) to get v_max and v_min, zero the spectrum between v_min and v_max,
  then produce two output cubes:
    * keep only velocities > global_vmax_global
    * keep only velocities < global_vmin_global

Usage example:
    python process_and_mask_cube.py input.fits
    python process_and_mask_cube.py input.fits -o out_pos.fits out_neg.fits

Notes:
- This version assumes spectral_cube is installed and usable.
- It expects a compiled rotation_curve_c module exposing calc_v_dev(...).
"""
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
import os
from spectral_cube import SpectralCube as sc
import time

# try import the C extension (or fallback to python implementation if you have it)
from rotation_model_numpy import calc_v_dev

def pix_to_galactic_l_b(celestial_wcs, xpix, ypix):
    lonlat = celestial_wcs.all_pix2world(xpix, ypix, 0)
    ctypes = celestial_wcs.wcs.ctype
    c0 = ctypes[0].upper()
    if "RA" in c0 or "DEC" in c0:
        sc = SkyCoord(ra=lonlat[0]*u.deg, dec=lonlat[1]*u.deg, frame='fk5')
        gal = sc.galactic
        return gal.l.deg, gal.b.deg
    else:
        return lonlat[0], lonlat[1]

def process_cube_mask_by_vdev(fits_in,
                             fits_out_pos_minvmax,
                             fits_out_neg_maxvmin,
                             h, r_gal, r_sun, v_sun, r_cut,
                             model,
                             v_dev):

    t0 = time.time()
    cube = sc.read(fits_in).with_spectral_unit(u.km / u.s)

    data = cube.unmasked_data[:].value  # (v, y, x)
    wcs = cube.wcs
    ny, nx = data.shape[-2], data.shape[-1]
    vel_axis = cube.spectral_axis.value

    pos_mask = np.zeros_like(data, dtype=bool)
    neg_mask = np.zeros_like(data, dtype=bool)

    # Create meshgrids for i (x) and j (y)
    i_grid, j_grid = np.meshgrid(np.arange(nx), np.arange(ny), indexing='ij')  # i_grid: (nx, ny), j_grid: (nx, ny)

    # Vectorized computation of Galactic coordinates
    l_deg, b_deg = pix_to_galactic_l_b(wcs.celestial, i_grid, j_grid)

    t1 = time.time()
    print(f"Loaded cube {fits_in}: shape={data.shape}, took {t1 - t0:.4f} s")

    # Vectorized call to calc_v_dev
    v_max, v_min = calc_v_dev(l_deg, b_deg,
                              h, r_gal, r_sun, v_sun, r_cut,
                              model=model, v_dev=v_dev)

    t2 = time.time()
    print(f"Computed v_max/v_min maps in {t2 - t1:.4f} s")

    vmax_map = v_max.T  # Transpose to (ny, nx)
    vmin_map = v_min.T

    # Create masks: broadcast vel_axis to (nv, ny, nx)
    pos_mask = vel_axis[:, None, None] > vmax_map[None, :, :]
    neg_mask = vel_axis[:, None, None] < vmin_map[None, :, :]

    t3 = time.time()
    print(f"Constructed masks in {t3 - t2:.4f} s")

    # 计算谱轴范围
    # 正速度: np.min(vmax_map) < v <= np.max(v)
    vmax_global = np.min(vmax_map)
    max_v = np.max(vel_axis)

    # 负速度: np.min(v) <= v < np.max(vmin_map)
    min_v = np.min(vel_axis)
    vmin_global = np.max(vmin_map)

    # 建立mask dask cube 以便保存（spectral_cube支持自带mask）
    pos_cube = sc(data=np.where(pos_mask, data, 0)*cube.unit, wcs=cube.wcs)
    neg_cube = sc(data=np.where(neg_mask, data, 0)*cube.unit, wcs=cube.wcs)

    pos_cube_slab = pos_cube.spectral_slab(vmax_global * u.km/u.s, max_v * u.km/u.s)
    neg_cube_slab = neg_cube.spectral_slab(min_v * u.km/u.s, vmin_global * u.km/u.s)

    t4 = time.time()
    print(f"Built masked cubes in {t4 - t3:.4f} s")

    # 直接用 spectral_cube 提供的 write 方法保存，header自动正确
    pos_cube_slab.write(fits_out_pos_minvmax, overwrite=True)
    neg_cube_slab.write(fits_out_neg_maxvmin, overwrite=True)

    t5 = time.time()
    print(f"Wrote output FITS files in {t5 - t4:.4f} s")

    print(f"Total time: {t5 - t0:.4f} s")

    return {
        "vmax_global": vmax_global,
        "vmin_global": vmin_global,
        "out_pos": fits_out_pos_minvmax,
        "out_neg": fits_out_neg_maxvmin
    }

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description="Mask FITS cube spectra: keep only > v_max or < v_min for each pixel; trim to global (min(vmax_map), max(v)) and (min(v), max(vmin_map)).")
    p.add_argument("fits_in", help="Input FITS cube")
    p.add_argument("-o", "--out", nargs=2, default=None, metavar=("OUT_POS","OUT_NEG"),
                   help="Two output FITS files: first keeps v>v_max, second keeps v<v_min. If not specified, defaults to <input>_+.fits and <input>_-.fits")
    p.add_argument("--height", type=float, default=5.0, help="Height above the Galactic plane in kpc (default 5)")
    p.add_argument("--r_gal", type=float, default=20.0, help="Maximum Galactocentric radius in kpc (default 20)")
    p.add_argument("--r_sun", type=float, default=8.5, help="Solar Galactocentric radius in kpc (default 8.5), only used in model=simple")
    p.add_argument("--v_sun", type=float, default=220.0, help="Solar orbital velocity in km/s (default 220), only used in model=simple")
    p.add_argument("--r_cut", type=float, default=0.5, help="Cutoff radius for solid body rotation in kpc (default 0.5), only used in model=simple")
    p.add_argument("--model", type=str, default="univ", choices=("simple","univ","linear","power"))
    p.add_argument("--vdev", type=float, default=50.0, help="dev passed into calc_v_dev (default 50)")
    args = p.parse_args()

    if args.out is None:
        base_name = os.path.splitext(args.fits_in)[0]
        out_pos = base_name + "_+.fits"
        out_neg = base_name + "_-.fits"
    else:
        out_pos, out_neg = args.out

    results = process_cube_mask_by_vdev(args.fits_in,
                                   out_pos,
                                   out_neg,
                                   h=args.height,
                                   r_gal=args.r_gal,
                                   r_sun=args.r_sun,
                                   v_sun=args.v_sun,
                                   r_cut=args.r_cut,
                                   model=args.model,
                                   v_dev=args.vdev)
    print("Done. Results:")
    for key, value in results.items():
        print(f"{key}: {value}")