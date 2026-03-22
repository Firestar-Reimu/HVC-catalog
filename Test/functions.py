# Data manipulation and analysis
import numpy as np
import pandas as pd

# Astropy modules
import astropy.units as u
from astropy.coordinates import SkyCoord

# Matplotlib for plotting
import matplotlib.pyplot as plt

# Scipy and related modules
import scipy.sparse as sp
from scipy.optimize import curve_fit

# Spectral cube handling
from spectral_cube import SpectralCube as sc

# Progress bar
from tqdm.notebook import tqdm

# Additional libraries
from regions import EllipseSkyRegion

# Multithreading
import os
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict

# Galaxy rotation curve
import sys

sys.path.append("../Galaxy Rotation Model/")

from rotation_model_numpy import calc_v_dev

# Suppress warnings from spectral_cube
import warnings
from spectral_cube.utils import PossiblySlowWarning

warnings.filterwarnings("ignore", category=PossiblySlowWarning)

# Set global options for numpy and matplotlib
np.seterr(divide="ignore", invalid="ignore")
np.set_printoptions(precision=10, suppress=True)
plt.rcParams.update({"figure.max_open_warning": 0, "font.size": 12})
delta_min = 1e-12


# Vectorized computation of Galactic coordinates (pixel to l,b)
def pix_to_galactic_l_b(celestial_wcs, xpix, ypix):
    lonlat = celestial_wcs.all_pix2world(xpix, ypix, 0)
    ctypes = celestial_wcs.wcs.ctype
    c0 = ctypes[0].upper()
    if "RA" in c0 or "DEC" in c0:
        sc = SkyCoord(ra=lonlat[0] * u.deg, dec=lonlat[1] * u.deg, frame="fk5")
        gal = sc.galactic
        return gal.l.deg, gal.b.deg
    else:
        return lonlat[0], lonlat[1]


# Precompute mapping cluster_id -> vyx rows (index rows)
def precompute_cluster_vyx(index: np.ndarray, labels: np.ndarray):
    """
    Precompute mapping cluster_id -> vyx rows (index rows).
    Returns a dict where keys are cluster IDs (labels >= 0) and values are arrays shape (m, 3).
    This does a single scan over `labels` and avoids repeated boolean indexing.
    """
    cluster_indices = defaultdict(list)
    for i, label in enumerate(labels):
        if label >= 0:
            cluster_indices[label].append(i)

    cluster_vyx = {}
    for cid, idx_list in cluster_indices.items():
        cluster_vyx[cid] = index[np.asarray(idx_list)]

    return cluster_vyx


# Calculate moment 0 using pure Python (faster for smaller clusters)
def moment_0_py(data, vyx, delta_v):
    # vyx is coordinates [v, y, x] of points in the cluster (labels == cluster_id)
    # Extract the velocity column and the coordinate column separately
    velo = vyx[:, 0]  # velocity [v]
    coords = vyx[:, 1:]  # coordinates [y, x]
    # Find unique rows in the x,y columns and their indices
    unique_coords, indices = np.unique(coords, axis=0, return_inverse=True)
    # unique_coords is a 2D array of unique [y, x] coordinates
    # indices is an array of the same length as coords, where each element is the index of the unique coordinate
    # This means that for each coordinate in coords, we can find several velocities in velo that correspond to it.

    # Calculate moment 0 of the velocity axis
    # for each unique pair of the coordinates
    def calc_moment_0(unique_coords, indices, velo, data):
        moment_0 = np.zeros(len(unique_coords))
        for i in range(len(unique_coords)):  # for each unique coordinate
            # Calculate the moment 0 for this coordinate
            # M0 = \int I_v dv = \sum I_v * delta_v
            m0 = 0
            coo = unique_coords[i]
            for v in velo[
                indices == i
            ]:  # sum over all velocities that correspond to this coordinate
                m0 += data[v, coo[0], coo[1]] * delta_v  # moment 0
            moment_0[i] = m0
        return moment_0

    moment_0 = calc_moment_0(unique_coords, indices, velo, data)
    # Combine the unique coordinates with their moment 0
    # wrap the result in a sparse matrix (COO format)
    # this is more memory efficient for large datasets
    moment_0_coo = sp.coo_array(
        (moment_0, (unique_coords[:, 0], unique_coords[:, 1])), shape=data[0].shape
    )
    return moment_0_coo


# Calculate moment 0 using pure Spectral-Cube (faster for larger clusters)
# https://github.com/radio-astro-tools/spectral-cube/blob/e98b6c3c05e3a21c6ca62524e1dea9582ad5cd38/spectral_cube/_moments.py#L170
def moment_0_sc(data, vyx, wcs):
    # create a cube which only contains the values of points in the cluster (labels == cluster_id)
    bool_array = np.zeros(data.shape, dtype=np.bool)
    bool_array[vyx[:, 0], vyx[:, 1], vyx[:, 2]] = (
        True  # mark all points in the cluster (labels == cluster_id) as True
    )
    cluster_data = (
        bool_array * data
    )  # multiply the data with the boolean array to get only the values of points in the cluster
    bool_cube = sc(cluster_data, wcs=wcs).with_spectral_unit(u.km / u.s)
    # Calculate moment 0 of the velocity axis
    moment_0 = bool_cube.moment(order=0).value
    # wrap the result in a sparse matrix (COO format)
    # this is more memory efficient for large datasets
    moment_0_coo = sp.coo_array(moment_0)
    return moment_0_coo


# Main function to calculate moment 0, using either pure Python or Spectral-Cube based on cluster size
def moment_0_func(data, cluster_vyx, n_clusters, delta_v, wcs):
    def compute_moment_0_for_cluster(cluster_id):
        vyx = cluster_vyx.get(cluster_id)
        size = len(vyx)
        if size < 100000:
            return moment_0_py(data, vyx, delta_v)
        else:
            return moment_0_sc(data, vyx, wcs)

    with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
        mapped = executor.map(compute_moment_0_for_cluster, range(n_clusters))
        moment_0_results = list(
            tqdm(
                mapped,
                total=n_clusters,
                desc="computing moment_0",
                mininterval=0.5,
                miniters=100,
            )
        )

    moment_0_cube = np.hstack(moment_0_results)
    return moment_0_cube


# Calculate moment 1 using pure Python (faster for smaller clusters)
def moment_1_py(data, vyx, velocities):
    # vyx is coordinates [v, y, x] of points in the cluster (labels == cluster_id)
    # Extract the velocity column and the coordinate column separately
    velo = vyx[:, 0]  # velocity [v]
    coords = vyx[:, 1:]  # coordinates [y, x]
    # Find unique rows in the x,y columns and their indices
    unique_coords, indices = np.unique(coords, axis=0, return_inverse=True)
    # unique_coords is a 2D array of unique [y, x] coordinates
    # indices is an array of the same length as coords, where each element is the index of the unique coordinate
    # This means that for each coordinate in coords, we can find several velocities in velo that correspond to it.

    # Calculate moment 1 of the velocity axis
    # for each unique pair of the coordinates
    def calc_moment_1(unique_coords, indices, velo, data):
        moment_1 = np.zeros(len(unique_coords))
        for i in range(len(unique_coords)):  # for each unique coordinate
            # Calculate the moment 1 for this coordinate
            # M1 = \int v I_v dv / M0 = \sum I_v * v / \sum I_v
            m0 = 0
            m1 = 0
            coo = unique_coords[i]
            for v in velo[
                indices == i
            ]:  # sum over all velocities that correspond to this coordinate
                m0 += data[
                    v, coo[0], coo[1]
                ]  # no need to multiply by delta_v since it cancels out in the division
                vel = velocities[v]
                m1 += data[v, coo[0], coo[1]] * vel  # moment 1
            moment_1[i] = m1 / m0
        return moment_1

    moment_1 = calc_moment_1(unique_coords, indices, velo, data)
    # Combine the unique coordinates with their moment 0
    # wrap the result in a sparse matrix (COO format)
    # this is more memory efficient for large datasets
    moment_1_coo = sp.coo_array(
        (moment_1, (unique_coords[:, 0], unique_coords[:, 1])), shape=data[0].shape
    )
    return moment_1_coo


# Calculate moment 0 using pure Spectral-Cube (faster for larger clusters)
# https://github.com/radio-astro-tools/spectral-cube/blob/e98b6c3c05e3a21c6ca62524e1dea9582ad5cd38/spectral_cube/_moments.py#L170
def moment_1_sc(data, vyx, wcs):
    # create a cube which only contains the values of points in the cluster (labels == cluster_id)
    bool_array = np.zeros(data.shape, dtype=np.bool)
    bool_array[vyx[:, 0], vyx[:, 1], vyx[:, 2]] = (
        True  # mark all points in the cluster (labels == cluster_id) as True
    )
    cluster_data = (
        bool_array * data
    )  # multiply the data with the boolean array to get only the values of points in the cluster
    bool_cube = sc(cluster_data, wcs=wcs).with_spectral_unit(u.km / u.s)
    # Calculate moment 1 of the velocity axis
    # Most pixels in the cube will be NaN due to the sparse nature of the data, which can cause RuntimeWarnings when calculating moments. We can suppress these warnings since we expect them and they don't indicate a problem with our data.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        moment_1 = bool_cube.moment(order=1).value
    # wrap the result in a sparse matrix (COO format)
    # this is more memory efficient for large datasets
    moment_1_coo = sp.coo_array(moment_1)
    return moment_1_coo


# Main function to calculate moment 1, using either pure Python or Spectral-Cube based on cluster size
def moment_1_func(data, cluster_vyx, n_clusters, velocities, wcs):
    def compute_moment_1_for_cluster(cluster_id):
        vyx = cluster_vyx.get(cluster_id)
        size = len(vyx)
        if size < 100000:
            return moment_1_py(data, vyx, velocities)
        else:
            return moment_1_sc(data, vyx, wcs)

    with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
        mapped = executor.map(compute_moment_1_for_cluster, range(n_clusters))
        moment_1_results = list(
            tqdm(
                mapped,
                total=n_clusters,
                desc="computing moment_1",
                mininterval=0.5,
                miniters=100,
            )
        )

    moment_1_cube = np.hstack(moment_1_results)
    return moment_1_cube


# Calculate moment 2 (FWHM) using pure Python (faster for smaller clusters)
def moment_2_py(data, vyx, velocities):
    # vyx is coordinates [v, y, x] of points in the cluster (labels == cluster_id)
    # Extract the velocity column and the coordinate column separately
    velo = vyx[:, 0]  # velocity [v]
    coords = vyx[:, 1:]  # coordinates [y, x]
    # Find unique rows in the x,y columns and their indices
    unique_coords, indices = np.unique(coords, axis=0, return_inverse=True)
    # unique_coords is a 2D array of unique [y, x] coordinates
    # indices is an array of the same length as coords, where each element is the index of the unique coordinate
    # This means that for each coordinate in coords, we can find several velocities in velo that correspond to it.

    # Calculate moment 2 of the velocity axis
    # for each unique pair of the coordinates
    def calc_moment_2(unique_coords, indices, velo, data):
        moment_2 = np.zeros(len(unique_coords))
        for i in range(len(unique_coords)):  # for each unique coordinate
            # Calculate the moment 1 for this coordinate
            # M2 = \int I_v (v - M1)^2 dv / M0 = \sum I_v * (v - M1)^2 / \sum I_v
            m0 = 0
            m1 = 0
            m2 = 0
            coo = unique_coords[i]
            # Calculate moment 0 and moment 1 for this coordinate
            for v in velo[
                indices == i
            ]:  # sum over all velocities that correspond to this coordinate
                m0 += data[
                    v, coo[0], coo[1]
                ]  # no need to multiply by delta_v since it cancels out in the division
                vel = velocities[v]
                m1 += data[v, coo[0], coo[1]] * vel  # moment 1
            m1 = m1 / m0
            # Calculate moment 2 for this coordinate
            for v in velo[indices == i]:
                vel = velocities[v]
                m2 += data[v, coo[0], coo[1]] * ((vel - m1) ** 2)  # moment 2
            m2 = m2 / m0
            # Convert moment 2 to FWHM
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("error")
                    sigma = np.sqrt(m2)  # linewidth_sigma
            except RuntimeWarning as e:
                warnings.simplefilter("ignore", category=RuntimeWarning)
                sigma = np.sqrt(m2)  # linewidth_sigma, warning (m2 < 0)
            fwhm = sigma * np.sqrt(8 * np.log(2))  # FWHM = sqrt(8 * log(2)) * sigma
            if np.abs(m2) > delta_min:
                moment_2[i] = fwhm
            else:
                moment_2[i] = 0
        moment_2 = np.nan_to_num(moment_2)
        return moment_2

    moment_2 = calc_moment_2(unique_coords, indices, velo, data)
    # Combine the unique coordinates with their moment 0
    # wrap the result in a sparse matrix (COO format)
    # this is more memory efficient for large datasets
    moment_2_coo = sp.coo_array(
        (moment_2, (unique_coords[:, 0], unique_coords[:, 1])), shape=data[0].shape
    )
    moment_2_coo.eliminate_zeros()  # eliminate zeros to save memory
    return moment_2_coo


# Calculate moment 2 using pure Spectral-Cube (faster for larger clusters)
# https://github.com/radio-astro-tools/spectral-cube/blob/e98b6c3c05e3a21c6ca62524e1dea9582ad5cd38/spectral_cube/_moments.py#L170
def moment_2_sc(data, vyx, wcs):
    # create a cube which only contains the values of points in the cluster (labels == cluster_id)
    bool_array = np.zeros(data.shape, dtype=np.bool)
    bool_array[vyx[:, 0], vyx[:, 1], vyx[:, 2]] = (
        True  # mark all points in the cluster (labels == cluster_id) as True
    )
    cluster_data = (
        bool_array * data
    )  # multiply the data with the boolean array to get only the values of points in the cluster
    bool_cube = sc(cluster_data, wcs=wcs).with_spectral_unit(u.km / u.s)
    # Calculate moment 2 (FWHM) of the velocity axis
    moment_2 = bool_cube.linewidth_fwhm().value
    moment_2 = np.nan_to_num(moment_2)
    moment_2[np.abs(moment_2) < delta_min] = 0  # set small values to 0
    # wrap the result in a sparse matrix (COO format)
    # this is more memory efficient for large datasets
    moment_2_coo = sp.coo_array(moment_2)
    return moment_2_coo


# Main function to calculate moment 2, using either pure Python or Spectral-Cube based on cluster size
def moment_2_func(data, cluster_vyx, n_clusters, velocities, wcs):
    def compute_moment_2_for_cluster(cluster_id):
        vyx = cluster_vyx.get(cluster_id)
        size = len(vyx)
        if size < 100000:
            return moment_2_py(data, vyx, velocities)
        else:
            return moment_2_sc(data, vyx, wcs)

    with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
        mapped = executor.map(compute_moment_2_for_cluster, range(n_clusters))
        moment_2_results = list(
            tqdm(
                mapped,
                total=n_clusters,
                desc="computing moment_2",
                mininterval=0.5,
                miniters=100,
            )
        )

    moment_2_cube = np.hstack(moment_2_results)
    return moment_2_cube


# Select clusters based on the criteria SNR, SIZE, and FWHM
def find_hvc_candidates(
    data,
    cluster_vyx,
    n_clusters,
    MEDIAN,
    STD,
    moment_0_cube,
    moment_2_cube,
    MIN_SNR,
    MIN_SIZE,
    MIN_FWHM,
):

    def compute_snr_size_fwhm_for_cluster(cluster_id):
        vyx = cluster_vyx[cluster_id]
        # SNR
        cluster_value = data[vyx[:, 0], vyx[:, 1], vyx[:, 2]]
        peak = np.max(cluster_value)
        snr = (peak - MEDIAN) / STD
        # SIZE
        moment_0_data = moment_0_cube[cluster_id].data
        size = len(np.argwhere(moment_0_data > 0))  # in pixels
        # FWHM
        moment_2_data = moment_2_cube[cluster_id].data
        if moment_2_data.size > 0:
            fwhm = np.mean(moment_2_data)
        else:
            # when all pixels in the cluster are in the same velocity channel,
            # moment_2_data.size == 0,
            # so the FWHM is set to 0 to avoid NaN values.
            fwhm = 0
        return np.array([snr, size, fwhm])

    # Parallel computation
    with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
        mapped = executor.map(compute_snr_size_fwhm_for_cluster, range(n_clusters))
        snr_size_fwhm_arr = list(
            tqdm(
                mapped,
                total=n_clusters,
                desc="computing snr_size_fwhm",
                mininterval=0.5,
                miniters=100,
            )
        )

    snr_size_fwhm_arr = np.array(snr_size_fwhm_arr)

    # Select HVC candidates based on criteria
    TF = (
        (snr_size_fwhm_arr[:, 0] > MIN_SNR)
        & (snr_size_fwhm_arr[:, 1] > MIN_SIZE)
        & (snr_size_fwhm_arr[:, 2] > MIN_FWHM)
    ).astype(int)

    hvc_candidates = np.nonzero(TF)[0]
    hvc_candidates_array = np.tile(TF, (3, 1)).T * snr_size_fwhm_arr

    SNR = hvc_candidates_array[hvc_candidates][:, 0]
    SIZE = hvc_candidates_array[hvc_candidates][:, 1]
    FWHM = hvc_candidates_array[hvc_candidates][:, 2]

    return hvc_candidates, SNR, SIZE, FWHM


# Calculate the minimum volume enclosing ellipse (MVEE) for a set of points, using the Khachiyan algorithm
def getMinVolEllipse(P, tolerance=1e-6, max_iter=10000):
    """
    Robust Khachiyan algorithm for MVEE.

    Parameters
    ----------
    P : ndarray, shape (N, d)
        Input points.
    tolerance : float
        Convergence tolerance (on max Mahalanobis violation).
    max_iter : int
        Safety cap on iterations.

    Returns
    -------
    xc, yc : float
        Ellipse center.
    a, b : float
        Semi-axis lengths (a >= b).
    theta : float
        Rotation angle in degrees (counter-clockwise, major axis).
    """

    P = np.asarray(P, dtype=float)
    N, d = P.shape

    # Build Q = [P^T; 1 ... 1]
    Q = np.vstack([P.T, np.ones(N)])

    # Initialize weights
    U = np.ones(N) / N

    for _ in range(max_iter):
        # V = Q diag(U) Q^T
        V = Q @ np.diag(U) @ Q.T

        # More stable than inv(V) @ Q
        V_inv_Q = np.linalg.solve(V, Q)
        M = Q.T @ V_inv_Q

        # Find most violating point
        diagM = np.diag(M)
        j = np.argmax(diagM)
        maximum = diagM[j]

        # --- Key improvement: correct stopping criterion ---
        if maximum <= d + 1 + tolerance:
            break

        # Khachiyan step size
        step = (maximum - d - 1) / ((d + 1) * (maximum - 1))
        step = np.clip(step, 0.0, 1.0)  # numerical safety

        U *= 1.0 - step
        U[j] += step

    # Center
    center = P.T @ U
    yc, xc = center  # note: swapped for image coordinates

    # Covariance-like matrix
    cov = P.T @ np.diag(U) @ P - np.outer(center, center)

    # Eigen-decomposition with numerical floor
    eigvals, eigvecs = np.linalg.eigh(cov)
    eigvals = np.maximum(eigvals, 1e-12)

    # A matrix of ellipse: (x-c)^T A (x-c) <= 1
    A = eigvecs @ np.diag(1.0 / (d * eigvals)) @ eigvecs.T

    # Extract geometric parameters
    lam, vec = np.linalg.eigh(A)
    idx = np.argsort(lam)  # ascending: smallest lambda = major axis
    lam = lam[idx]
    vec = vec[:, idx]

    a = 1.0 / np.sqrt(lam[0])  # major semi-axis
    b = 1.0 / np.sqrt(lam[1])  # minor semi-axis

    # Angle of major axis
    theta = np.degrees(np.arctan2(vec[1, 0], vec[0, 0]))
    theta = theta + 90  # rotate to image coordinates

    return xc, yc, a, b, theta


# Calculate the FWHM ellipse parameters (center, semi-major axis, semi-minor axis, angle) for a given moment 0 crop
def fwhm_ellipse(moment_0_crop):
    points_above_half_max = np.argwhere(moment_0_crop > 0.5 * np.max(moment_0_crop))
    points_above_zero = np.argwhere(moment_0_crop > 0)
    area = len(points_above_half_max)  # Number of points above half maximum
    # Calculate the minimum volume ellipse for the included points
    try:
        # First attempt calculation using points_above_half_max
        xc, yc, a, b, theta = getMinVolEllipse(points_above_half_max)
    except np.linalg.LinAlgError as e:
        # If the above calculation fails (e.g., insufficient points, singular matrix, etc.), catch the exception
        print(f"Singular Matrix Error: {e}. Using all points for ellipse fitting.")
        # Use the alternative points dataset for calculation
        xc, yc, a, b, theta = getMinVolEllipse(points_above_zero)
    return xc, yc, a, b, theta, area


# Calculate the size of the cluster based on the moment 0 cube
def calc_size(cluster_id, moment_0_cube, wcs):
    # Get the moment 0 array from the sparse matrix
    moment_0 = moment_0_cube[cluster_id].toarray()
    # Find the non-zero indices in the moment 0 cube
    nonzero_indices = np.nonzero(moment_0)
    y_min, y_max = np.min(nonzero_indices[0]), np.max(nonzero_indices[0])
    x_min, x_max = np.min(nonzero_indices[1]), np.max(nonzero_indices[1])
    moment_0_crop = moment_0[y_min : y_max + 1, x_min : x_max + 1]
    height, width = moment_0_crop.shape
    # Calculate the extent of the cropped moment 0 array
    # note that RA is left bigger, right smaller, and DEC is bottom smaller, top bigger
    extent = [
        wcs.celestial.pixel_to_world_values(x_min, 0)[0],
        wcs.celestial.pixel_to_world_values(x_max, 0)[0],
        wcs.celestial.pixel_to_world_values(0, y_min)[1],
        wcs.celestial.pixel_to_world_values(0, y_max)[1],
    ]
    x_min, x_max, y_min, y_max = extent  # now extent is in world coordinates
    xc, yc, a, b, theta, area = fwhm_ellipse(moment_0_crop)
    x_scale = (x_max - x_min) / (width - 1)  # pixel-to-world scale factor, here 0.025
    y_scale = (y_max - y_min) / (height - 1)  # pixel-to-world scale factor, here 0.025
    xc = x_min + (xc * x_scale)
    yc = y_min + (yc * y_scale)
    area = area * np.abs(x_scale * y_scale)  # area from pixel**2 to deg**2
    a = a * np.abs(x_scale)  # a is the semi-major axis in degrees
    b = b * np.abs(y_scale)  # b is the semi-minor axis in degrees
    return moment_0_crop, extent, xc, yc, a, b, theta, area


# Find the peak between v_min and v_max
def find_peak_in_range(spectrum, velocity_axis, v_min, v_max):
    # Find the index range corresponding to v_min and v_max
    idx_range_start = np.argmin(np.abs(velocity_axis - v_min))
    idx_range_end = np.argmin(np.abs(velocity_axis - v_max)) + 1

    # Find the peak within the specified range
    spectrum_in_range = spectrum[idx_range_start:idx_range_end]
    idx_peak_local = np.argmax(spectrum_in_range)
    idx_peak = idx_range_start + idx_peak_local
    v_peak = velocity_axis[idx_peak]

    return idx_peak, v_peak


# Extract CRAFTS and HI4PI spectra for HVC candidates using the same elliptical sky region.
def extract_candidate_spectra(
    hvc_candidates,
    *,
    labels,
    index,
    data,
    velocities,
    cube,
    ellipse_params,
    hi4pi_cube,
    hi4pi_velocities,
    wcs,
    find_peak_in_range,
    extract_spectra_around_peak,
    half_width=50.0,
):
    """
    Extract CRAFTS and HI4PI spectra for HVC candidates using the same elliptical sky region.

    Parameters
    ----------
    hvc_candidates : array-like
        Candidate cluster IDs (values that appear in `labels`).
    labels : ndarray
        Label cube, same shape as `data`, where labels == cluster_id indicates membership.
    index : ndarray, shape (Npoints, 3)
        Array of voxel coordinates [v, y, x] corresponding to flattened/selected points.
        Used as: vyx = index[labels == cluster_id]
    data : ndarray
        Data cube (e.g., brightness temperature), shape (nv, ny, nx).
    velocities : ndarray
        Velocity axis for `cube` / `data`, indexed by v channel.
    cube : SpectralCube-like
        CRAFTS spectral cube; must support subcube_from_regions and slicing cube[:, y, x].value
    ellipse_params : ndarray, shape (len(hvc_candidates), 5)
        For each candidate i: (xc, yc, a, b, theta)
        where xc,yc in deg FK5; a,b in deg; theta in deg.
    hi4pi_cube : SpectralCube-like
        HI4PI spectral cube; must support subcube_from_regions and slicing hi4pi_cube[:, y, x].value
    hi4pi_velocities : ndarray
        Velocity axis for HI4PI.
    wcs : astropy.wcs.WCS
        WCS for the CRAFTS cube (used to convert peak pixel->world when HI4PI fallback is needed).
    find_peak_in_range : callable
        Function: (spectrum, velocities, v_min, v_max) -> (idx_peak, v_peak)
    extract_spectra_around_peak : callable
        Function:
          (crafts_spectrum, hi4pi_spectrum, velocities, hi4pi_velocities, v_min, v_max, v_peak, half_width)
          -> (crafts_cut, crafts_vel_cut, hi4pi_cut, hi4pi_vel_cut)
    half_width : float
        Half width (km/s) around v_peak for cutting spectra.

    Returns
    -------
    crafts_spectra, hi4pi_spectra : np.ndarray (n_candidates, n_chan_cut)
    crafts_velocity_axis_list, hi4pi_velocity_axis_list : np.ndarray (n_candidates, n_chan_cut)
    VPEAK, VMIN, VMAX : np.ndarray (n_candidates,)
    """
    crafts_spectra = []
    hi4pi_spectra = []
    crafts_velocity_axis_list = []
    hi4pi_velocity_axis_list = []
    VPEAK = []
    VMIN = []
    VMAX = []

    for i, cluster_id in enumerate(hvc_candidates):
        # Coordinates [v, y, x] of points in this cluster
        vyx = index[labels == cluster_id]

        # If a candidate has no voxels (shouldn't happen, but safe-guard)
        if vyx.size == 0:
            raise ValueError(f"Empty cluster: i={i}, cluster_id={cluster_id}")

        # Global v_min and v_max for the cluster (by channel indices)
        v_min = velocities[vyx[:, 0].min()]
        v_max = velocities[vyx[:, 0].max()]

        # Build elliptical region (FK5 deg)
        xc, yc, a, b, theta = ellipse_params[i, :]
        center = SkyCoord(xc, yc, unit="deg", frame="fk5")
        region = EllipseSkyRegion(
            center=center,
            width=a * u.deg,
            height=b * u.deg,
            angle=-theta * u.deg,
        )

        # --- Extract CRAFTS spectrum ---
        crafts_subcube = cube.subcube_from_regions([region])
        crafts_subcube_data = crafts_subcube.unmasked_data[:, :, :].value

        peak_vyx = None  # cache for fallbacks
        if crafts_subcube_data.size == 0:
            # Fallback: use spectrum at peak voxel of this cluster (in CRAFTS pixel space)
            cluster_value = data[vyx[:, 0], vyx[:, 1], vyx[:, 2]]
            peak_idx = int(np.argmax(cluster_value))
            peak_vyx = vyx[peak_idx]  # [v, y, x]
            y_peak, x_peak = int(peak_vyx[1]), int(peak_vyx[2])
            crafts_spectrum = cube[:, y_peak, x_peak].value
        else:
            crafts_spectrum = np.mean(crafts_subcube_data, axis=(1, 2))

        # --- Extract HI4PI spectrum ---
        hi4pi_subcube = hi4pi_cube.subcube_from_regions([region])
        hi4pi_subcube_data = hi4pi_subcube.unmasked_data[:, :, :].value

        if hi4pi_subcube_data.size == 0:
            # Fallback: map CRAFTS peak position to world -> HI4PI pixel
            if peak_vyx is None:
                cluster_value = data[vyx[:, 0], vyx[:, 1], vyx[:, 2]]
                peak_idx = int(np.argmax(cluster_value))
                peak_vyx = vyx[peak_idx]

            # pixel_to_world_values expects (x, y, v) or (x, y, z) depending on WCS;
            # keep consistent with your original usage.
            peak_world = wcs.pixel_to_world_values(
                float(peak_vyx[2]), float(peak_vyx[1]), float(peak_vyx[0])
            )
            peak_sky = SkyCoord(peak_world[0], peak_world[1], unit="deg", frame="fk5")

            hi4pi_wcs = hi4pi_cube.wcs.celestial
            x_hi4pi, y_hi4pi = hi4pi_wcs.world_to_pixel(peak_sky)
            hi4pi_spectrum = hi4pi_cube[:, int(y_hi4pi), int(x_hi4pi)].value
        else:
            hi4pi_spectrum = np.mean(hi4pi_subcube_data, axis=(1, 2))

        # --- Peak and cut around peak ---
        _, v_peak = find_peak_in_range(crafts_spectrum, velocities, v_min, v_max)

        crafts_cut, crafts_vel_cut, hi4pi_cut, hi4pi_vel_cut = (
            extract_spectra_around_peak(
                crafts_spectrum,
                hi4pi_spectrum,
                velocities,
                hi4pi_velocities,
                v_min,
                v_max,
                v_peak,
                half_width=half_width,
            )
        )

        crafts_spectra.append(crafts_cut)
        hi4pi_spectra.append(hi4pi_cut)
        crafts_velocity_axis_list.append(crafts_vel_cut)
        hi4pi_velocity_axis_list.append(hi4pi_vel_cut)
        VPEAK.append(v_peak)
        VMIN.append(v_min)
        VMAX.append(v_max)

    return (
        np.array(crafts_spectra, dtype=object),
        np.array(hi4pi_spectra, dtype=object),
        np.array(crafts_velocity_axis_list, dtype=object),
        np.array(hi4pi_velocity_axis_list, dtype=object),
        np.array(VPEAK, dtype=float),
        np.array(VMIN, dtype=float),
        np.array(VMAX, dtype=float),
    )


# Extract spectra around the peak velocity
# half width is 50 km/s by default
# if v_min or v_max is outside half width, use v_min or v_max instead
def extract_spectra_around_peak(
    crafts_spectrum,
    hi4pi_spectrum,
    velocities,
    hi4pi_velocities,
    v_min,
    v_max,
    peak,
    half_width=50.0,
):
    # Determine final velocity window
    v_left = min(peak - half_width, v_min)
    v_right = max(peak + half_width, v_max)

    # CRAFTS
    mask_crafts = (velocities >= v_left) & (velocities <= v_right)
    crafts_cut = crafts_spectrum[mask_crafts]
    crafts_vel_cut = velocities[mask_crafts]

    # HI4PI
    mask_hi4pi = (hi4pi_velocities >= v_left) & (hi4pi_velocities <= v_right)
    hi4pi_cut = hi4pi_spectrum[mask_hi4pi]
    hi4pi_vel_cut = hi4pi_velocities[mask_hi4pi]

    return crafts_cut, crafts_vel_cut, hi4pi_cut, hi4pi_vel_cut


def fit_gaussian(velocities, spectrum, v_peak, window=15.0):
    def gaussian(v, amp, v0, sigma):
        return amp * np.exp(-0.5 * ((v - v0) / sigma) ** 2)

    # 1. Fitting window
    mask = (velocities >= v_peak - window) & (velocities <= v_peak + window)
    v_fit = velocities[mask]
    y_fit = spectrum[mask]

    if len(v_fit) < 5:
        raise RuntimeError("Too few points in fitting window.")

    # 2. Initial parameter estimates
    amp0 = y_fit.max()
    v0_0 = v_fit[y_fit.argmax()]
    sigma0 = 8.5  # initial guess for sigma, corresponds to FWHM ~ 20 km/s

    p0 = (amp0, v0_0, sigma0)

    # 3. Gaussian Fitting
    popt, pcov = curve_fit(gaussian, v_fit, y_fit, p0=p0, maxfev=100000)
    amp, v_peak_fit, sigma = popt
    fwhm = 2 * np.sqrt(2 * np.log(2)) * sigma

    # 4. Calculate RMS of residuals
    y_model = gaussian(v_fit, *popt)
    residuals = y_fit - y_model
    rms_residual = (residuals**2).mean() ** 0.5

    # print(f"amp={amp:.3f} K, v_peak={v_peak_fit:.3f} km/s, FWHM={fwhm:.3f} km/s")
    return amp, v_peak_fit, fwhm, v_fit, y_model, rms_residual


# RA to HMS, DEC to DMS
def ra_dec_to_hms_dms(ra, dec):
    c = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="fk5")
    ra_hms = c.ra.to_string(unit=u.hour, sep=":", precision=0, pad=True)
    dec_dms = c.dec.to_string(
        unit=u.degree, sep=":", precision=0, pad=True, alwayssign=True
    )
    return ra_hms, dec_dms


# Build HVC catalog DataFrame
def build_hvc_catalog(
    hvc_candidates,
    *,
    ID,
    ellipse_params,
    VPEAK,
    SNR,
    FWHM,
    MEDIAN,
    STD,
    moment_0_cube,
    ra_dec_to_hms_dms,
    calc_v_dev,
    vdev_model="poly",
):
    cols = [
        "ID",
        "RA",
        "DEC",
        "GLON",
        "GLAT",
        "VLSR",
        "VGSR",
        "VDEV",
        "SIZEa",
        "SIZEb",
        "SNR",
        "FWHM",
        "TPEAK",
        "N_HI",
    ]
    catalog = pd.DataFrame(columns=cols)

    for i in range(len(hvc_candidates)):
        cluster_id = hvc_candidates[i]

        # HVC ID
        hvc_id = ID[i]

        # ellipse parameters
        xc, yc, a, b, theta = ellipse_params[i]
        a_arcmin = 2 * a * 60  # major axis diameter in arcmin
        b_arcmin = 2 * b * 60  # minor axis diameter in arcmin

        # coordinates
        c_eq = SkyCoord(xc, yc, frame="fk5", unit="deg")
        c_gal = c_eq.galactic
        ra = c_eq.ra.value
        dec = c_eq.dec.value
        ra_hms, dec_dms = ra_dec_to_hms_dms(ra, dec)
        glon = c_gal.l.value
        glat = c_gal.b.value

        # peak brightness temperature derived from SNR
        snr = SNR[i]
        t_peak = snr * STD + MEDIAN

        # velocities
        vlsr = VPEAK[i]
        vmodel = calc_v_dev(glon, glat, model=vdev_model)
        vdev = vlsr - vmodel[0]
        vgsr = vlsr + 220 * np.sin(np.deg2rad(glon)) * np.cos(np.deg2rad(glat))

        # HI column density (1e20 cm^-2)
        n_hi = np.max(moment_0_cube[cluster_id]) * 1.823 / 100

        # FWHM
        fwhm = FWHM[i]

        catalog.loc[i] = {
            "ID": hvc_id,
            "RA": ra_hms,
            "DEC": dec_dms,
            "GLON": glon,
            "GLAT": glat,
            "VLSR": vlsr,
            "VGSR": vgsr,
            "VDEV": vdev,
            "SIZEa": a_arcmin,
            "SIZEb": b_arcmin,
            "SNR": snr,
            "FWHM": fwhm,
            "TPEAK": t_peak,
            "N_HI": n_hi,
        }
    return catalog
