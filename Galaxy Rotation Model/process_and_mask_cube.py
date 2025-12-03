#!/usr/bin/env python3
"""
process_and_mask_cube.py

Simplified processing script:

- Uses spectral_cube exclusively to extract the spectral axis in km/s:
    from spectral_cube import SpectralCube as sc
    cube = sc.read(file).with_spectral_unit(u.km / u.s)
    velocities = cube.spectral_axis.value

- For each spatial pixel (x,y) compute Galactic coordinates (l,b) using the cube WCS,
  call calc_v_dev(l,b,...) to get v_max and v_min, zero the spectrum between v_min..v_max,
  then produce two output cubes:
    * keep only velocities > global_min_vmax
    * keep only velocities < global_max_vmin

Usage example:
    python process_and_mask_cube.py input.fits -o out_pos.fits out_neg.fits

Notes:
- This version assumes spectral_cube is installed and usable.
- It expects a compiled rotation_curve_c module exposing calc_v_dev(...).
"""
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
import numpy as np
from tqdm import tqdm
import warnings
from astropy.coordinates import SkyCoord

# Use spectral_cube as requested
from spectral_cube import SpectralCube as sc

# try import the C extension (or fallback to python implementation if you have it)
try:
    #import rotation_curve_c as rcmod
    #calc_v_dev = rcmod.calc_v_dev
    from Galaxy_Rotation_Model import calc_v_dev_py as calc_v_dev
except Exception:
    raise ImportError("rotation_curve_c module not found. Build and install the C extension or provide a Python calc_v_dev and modify the script.")

def load_cube(path, hdu=0, memmap=True):
    hdul = fits.open(path, memmap=memmap)
    hduobj = hdul[hdu]
    header = hduobj.header
    data = hduobj.data
    if data is not None:
        data = np.asarray(data)
    try:
        wcs = WCS(header)
        if getattr(wcs.wcs, "naxis", 0) == 0:
            wcs = None
    except Exception:
        wcs = None
    hdul.close()
    return data, header, wcs

def get_spectral_velocity_axis(fits_path, nchan, rest_freq=None):
    """
    Simplified: use spectral_cube to read the file and convert spectral axis to km/s.
    Returns a numpy array of length nchan with velocities in km/s.

    This intentionally implements only the single path you specified.
    """
    cube = sc.read(fits_path)
    # Convert to km/s (spectral_cube will try to interpret the axis)
    if rest_freq is not None:
        rest_q = u.Quantity(rest_freq, u.Hz)
        cube_kms = cube.with_spectral_unit(u.km / u.s, rest_value=rest_q)
    else:
        cube_kms = cube.with_spectral_unit(u.km / u.s)
    vel_axis = cube_kms.spectral_axis.value
    vel_axis = np.asarray(vel_axis, dtype=float)
    # If lengths mismatch, warn but still return what spectral_cube provided
    if len(vel_axis) != nchan:
        warnings.warn(f"spectral_cube returned {len(vel_axis)} channels but expected {nchan}.")
    return vel_axis

def pix_to_galactic_l_b(celestial_wcs, xpix, ypix):
    """
    Given a 2D celestial_wcs and pixel indices (x,y) in numpy indexing
    (i.e., x is column index, y is row index), return galactic (l,b) in degrees.
    Uses FK5 for RA/Dec interpretation (per request).
    """
    if celestial_wcs is None:
        return None, None
    try:
        lonlat = celestial_wcs.all_pix2world(float(xpix), float(ypix), 0)
        ctypes = getattr(celestial_wcs.wcs, "ctype", [])
        if ctypes:
            c0 = ctypes[0].upper()
            # If header gives RA/Dec, convert using FK5 (requested) to Galactic
            if "RA" in c0 or "DEC" in c0 or ("GLON" not in c0 and "GLAT" not in c0):
                ra = lonlat[0]
                dec = lonlat[1]
                sc = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='fk5')
                gal = sc.galactic
                return float(gal.l.deg), float(gal.b.deg)
            else:
                return float(lonlat[0]), float(lonlat[1])
        else:
            return float(lonlat[0]), float(lonlat[1])
    except Exception:
        return None, None

def find_spectral_index(wcs):
    ctypes = getattr(wcs.wcs, "ctype", [])
    for i, ct in enumerate(ctypes):
        up = ct.upper()
        if any(sub in up for sub in ("FREQ", "VELO", "VRAD", "VOPT", "WAVE", "AWAV", "BETA")):
            return i
    return None

def process_cube_mask_by_vdev(fits_in,
                             fits_out_gt_minvmax,
                             fits_out_lt_maxvmin,
                             hdu=0,
                             spec_axis=0,
                             model='univ',
                             v_dev=50.0,
                             rest_freq=None,
                             memmap=True):
    """
    Main processing function using simplified spectral axis extraction (spectral_cube).
    """
    if rest_freq is not None:
        rest_freq = u.Quantity(rest_freq, u.Hz)

    data, header, wcs = load_cube(fits_in, hdu=hdu, memmap=memmap)
    if data is None:
        raise ValueError("Input FITS has no data.")

    # Move spectral axis to axis 0
    if spec_axis != 0:
        data = np.moveaxis(data, spec_axis, 0)
    if data.ndim < 3:
        raise ValueError("Data must be a cube with at least 3 axes (spec,y,x).")
    nz = data.shape[0]
    ny = data.shape[-2]
    nx = data.shape[-1]

    # get spectral velocity axis (1D array length nz) in km/s using spectral_cube
    vel_axis = get_spectral_velocity_axis(fits_in, nz, rest_freq=rest_freq.to(u.Hz).value)
    if vel_axis is None:
        raise RuntimeError("Failed to determine velocity axis via spectral_cube.")

    masked_cube = data.copy()
    vmax_map = np.full((ny, nx), np.nan, dtype=float)
    vmin_map = np.full((ny, nx), np.nan, dtype=float)

    # Prepare celestial WCS for pixel->sky (spatial) if possible
    celestial_wcs = None
    if wcs is not None:
        try:
            celestial_wcs = wcs.celestial
        except Exception:
            celestial_wcs = None

    total = nx * ny
    pbar = tqdm(total=total, desc="Masking pixels")
    for j in range(ny):  # y (row)
        for i in range(nx):  # x (col)
            if celestial_wcs is not None:
                l_deg, b_deg = pix_to_galactic_l_b(celestial_wcs, i, j)
            else:
                l_deg, b_deg = None, None

            if l_deg is None or b_deg is None:
                vmax_map[j, i] = np.nan
                vmin_map[j, i] = np.nan
            else:
                # Attempt calling calc_v_dev using the Python-default signature (no n_sample):
                # calc_v_dev(l, b, h=5, r_gal=20, r_sun=8.5, v_sun=220, model='univ', v_dev=50)
                try:
                    v_max, v_min = calc_v_dev(float(l_deg), float(b_deg),
                                               h=5, r_gal=20, r_sun=8.5, v_sun=220,
                                               model=model, v_dev=v_dev)
                except Exception as e:
                    warnings.warn(f"calc_v_dev raised at pixel ({i},{j}) l={l_deg},b={b_deg}: {e}")
                    v_max, v_min = np.nan, np.nan

                vmax_map[j, i] = float(v_max) if np.isfinite(v_max) else np.nan
                vmin_map[j, i] = float(v_min) if np.isfinite(v_min) else np.nan

                if np.isfinite(v_min) and np.isfinite(v_max):
                    vlow = min(v_min, v_max)
                    vhigh = max(v_min, v_max)
                    mask = (vel_axis >= vlow) & (vel_axis <= vhigh)
                    if mask.any():
                        masked_cube[mask, j, i] = 0.0
            pbar.update(1)
    pbar.close()

    if np.all(np.isnan(vmax_map)):
        raise RuntimeError("All vmax values are NaN: cannot compute global thresholds.")
    global_min_vmax = np.nanmin(vmax_map)
    global_max_vmin = np.nanmax(vmin_map)

    cube_gt_minvmax = masked_cube.copy()
    cube_lt_maxvmin = masked_cube.copy()

    mask_le = vel_axis <= global_min_vmax
    if mask_le.any():
        cube_gt_minvmax[mask_le, :, :] = 0.0

    mask_ge = vel_axis >= global_max_vmin
    if mask_ge.any():
        cube_lt_maxvmin[mask_ge, :, :] = 0.0

    out_header = header.copy()
    if spec_axis != 0:
        cube_to_write_gt = np.moveaxis(cube_gt_minvmax, 0, spec_axis)
        cube_to_write_lt = np.moveaxis(cube_lt_maxvmin, 0, spec_axis)
    else:
        cube_to_write_gt = cube_gt_minvmax
        cube_to_write_lt = cube_lt_maxvmin

    fits.PrimaryHDU(data=cube_to_write_gt, header=out_header).writeto(fits_out_gt_minvmax, overwrite=True)
    fits.PrimaryHDU(data=cube_to_write_lt, header=out_header).writeto(fits_out_lt_maxvmin, overwrite=True)

    return {
        "vmax_map": vmax_map,
        "vmin_map": vmin_map,
        "global_min_vmax": global_min_vmax,
        "global_max_vmin": global_max_vmin,
        "out_gt": fits_out_gt_minvmax,
        "out_lt": fits_out_lt_maxvmin
    }

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description="Mask FITS cube spectra between v_min and v_max computed per-pixel from Galactic (l,b).")
    p.add_argument("fits_in", help="Input FITS cube")
    p.add_argument("-o", "--out", nargs=2, required=True, metavar=("OUT_POS","OUT_NEG"),
                   help="Two output FITS files: first keeps velocities > min(v_max), second keeps velocities < max(v_min)")
    p.add_argument("--hdu", type=int, default=0)
    p.add_argument("--spec-axis", type=int, default=0, help="numpy axis index that is spectral (default 0 => data.shape=(spec,y,x))")
    p.add_argument("--model", type=str, default="univ", choices=("simple","univ","linear","linear2","flat","power"))
    p.add_argument("--vdev", type=float, default=50.0, help="dev passed into calc_v_dev (default 0)")
    p.add_argument("--rest-freq", type=float, default=1.420405751770e9, help="rest frequency in Hz (optional, for freq->velocity conversion)")
    args = p.parse_args()

    out_pos, out_neg = args.out
    res = process_cube_mask_by_vdev(args.fits_in,
                                   out_pos,
                                   out_neg,
                                   hdu=args.hdu,
                                   spec_axis=args.spec_axis,
                                   model=args.model,
                                   v_dev=args.vdev,
                                   rest_freq=u.Quantity(args.rest_freq, u.Hz),
                                   memmap=True)
    print("Done. Results:", res)