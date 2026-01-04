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
    python process_and_mask_cube.py ./xx/yy/*.fits

Notes:
- This version assumes spectral_cube is installed and usable.
- It expects a compiled rotation_curve_c module exposing calc_v_dev(...).
"""
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
import os
from spectral_cube import SpectralCube as sc

# import calc_v_dev from another python file
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

    cube = sc.read(fits_in).with_spectral_unit(u.km / u.s)

    data = cube.unmasked_data[:].value  # (v, y, x)
    wcs = cube.wcs
    ny, nx = data.shape[-2], data.shape[-1]
    vel_axis = cube.spectral_axis.value

    pos_mask = np.zeros_like(data, dtype=bool)
    neg_mask = np.zeros_like(data, dtype=bool)

    # Create meshgrids for i (x) and j (y)
    i_grid, j_grid = np.meshgrid(np.arange(nx), np.arange(ny), indexing='ij')

    # Vectorized computation of Galactic coordinates
    l_deg, b_deg = pix_to_galactic_l_b(wcs.celestial, i_grid, j_grid)

    # Vectorized call to calc_v_dev
    v_max, v_min = calc_v_dev(l_deg, b_deg,
                              h, r_gal, r_sun, v_sun, r_cut,
                              model=model, v_dev=v_dev)

    vmax_map = v_max.T  # Transpose to (ny, nx)
    vmin_map = v_min.T

    # Create masks: broadcast vel_axis to (nv, ny, nx)
    pos_mask = vel_axis[:, None, None] > vmax_map[None, :, :]
    neg_mask = vel_axis[:, None, None] < vmin_map[None, :, :]

    # Calculate global vmin and vmax for trimming
    # Positive velocities: np.min(vmax_map) < v <= np.max(v)
    vmax_global = np.min(vmax_map)
    max_v = np.max(vel_axis)

    # Negative velocities: np.min(v) <= v < np.max(vmin_map)
    min_v = np.min(vel_axis)
    vmin_global = np.max(vmin_map)

    # Build masked dask cubes for saving (spectral_cube supports built-in masks)
    pos_cube = sc(data=np.where(pos_mask, data, 0)*cube.unit, wcs=cube.wcs)
    neg_cube = sc(data=np.where(neg_mask, data, 0)*cube.unit, wcs=cube.wcs)

    pos_cube_slab = pos_cube.spectral_slab(vmax_global * u.km/u.s, max_v * u.km/u.s)
    neg_cube_slab = neg_cube.spectral_slab(min_v * u.km/u.s, vmin_global * u.km/u.s)

    # Use spectral_cube's built-in write method to save, header is automatically correct
    pos_cube_slab.write(fits_out_pos_minvmax, overwrite=True)
    neg_cube_slab.write(fits_out_neg_maxvmin, overwrite=True)

    return {
        "vmax_global": vmax_global,
        "vmin_global": vmin_global,
        "out_pos": fits_out_pos_minvmax,
        "out_neg": fits_out_neg_maxvmin
    }

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description="Mask FITS cube spectra: keep only > v_max or < v_min for each pixel; trim to global (min(vmax_map), max(v)) and (min(v), max(vmin_map)).")
    p.add_argument("fits_in", nargs='+', help="Input FITS cube(s)")
    p.add_argument("-o", "--out", nargs=2, default=None, metavar=("OUT_POS","OUT_NEG"),
                   help="Two output FITS files: first keeps v>v_max, second keeps v<v_min. If not specified, defaults to <input>_+.fits and <input>_-.fits. For multiple inputs, outputs are generated per input.")
    p.add_argument("--height", type=float, default=5.0, help="Height above the Galactic plane in kpc (default 5)")
    p.add_argument("--r_gal", type=float, default=20.0, help="Maximum Galactocentric radius in kpc (default 20)")
    p.add_argument("--r_sun", type=float, default=8.5, help="Solar Galactocentric radius in kpc (default 8.5), only used in model=simple")
    p.add_argument("--v_sun", type=float, default=220.0, help="Solar orbital velocity in km/s (default 220), only used in model=simple")
    p.add_argument("--r_cut", type=float, default=0.5, help="Cutoff radius for solid body rotation in kpc (default 0.5), only used in model=simple")
    p.add_argument("--model", type=str, default="simple", choices=("simple","univ","linear","poly"))
    p.add_argument("--v_dev", type=float, default=0.0, help="dev passed into calc_v_dev (default 50)")
    args = p.parse_args()

    for fits_in in args.fits_in:
        if args.out is None:
            base_name = os.path.splitext(fits_in)[0]
            out_pos = base_name + "_+.fits"
            out_neg = base_name + "_-.fits"
        else:
            # For multiple files, if -o is specified, it might not make sense; perhaps warn or skip
            if len(args.fits_in) > 1:
                print("Warning: -o specified for multiple input files. Using default naming for each.")
                base_name = os.path.splitext(fits_in)[0]
                out_pos = base_name + "_+.fits"
                out_neg = base_name + "_-.fits"
            else:
                out_pos, out_neg = args.out

        # print current arguments
        print(f"Arguments for {fits_in}: height: {args.height}, r_gal: {args.r_gal}, r_sun: {args.r_sun}, v_sun: {args.v_sun}, r_cut: {args.r_cut}, model: {args.model}, v_dev: {args.v_dev}")

        results = process_cube_mask_by_vdev(fits_in,
                                       out_pos,
                                       out_neg,
                                       h=args.height,
                                       r_gal=args.r_gal,
                                       r_sun=args.r_sun,
                                       v_sun=args.v_sun,
                                       r_cut=args.r_cut,
                                       model=args.model,
                                       v_dev=args.v_dev)
        print(f"Processed {fits_in}. Results:")
        for key, value in results.items():
            print(f"{key}: {value}")
        print("-" * 50)