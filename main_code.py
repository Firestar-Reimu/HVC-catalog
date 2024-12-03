# %%
import numpy as np
import scipy.sparse as sp
from scipy.ndimage import gaussian_filter
from scipy.signal import find_peaks, peak_prominences
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import fits
from astropy.stats import sigma_clip
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.visualization import simple_norm
from astropy.wcs import WCS
from spectral_cube import SpectralCube as sc
from sklearn.cluster import DBSCAN
from tqdm import trange
import os
import warnings
from copy import deepcopy

np.seterr(divide="ignore", invalid="ignore")
np.set_printoptions(precision=10, suppress=True)
plt.rcParams.update({'figure.max_open_warning': 0})
delta_min = 1e-12

# %matplotlib widget

# %%
# 从data cube导入数据

file = "HVC_229_234_4_14.fits"
hdul = fits.open(file)
hdr = hdul[0].header
wcs = WCS(hdr)
data_nan = hdul[0].data.astype(np.float64)
pixel_to_arcmin = (hdr["CDELT2"] * u.deg).to(u.arcmin).value  # 1 pixel = 1.5 arcmin

MEAN = np.nanmean(data_nan)
STD = np.nanstd(data_nan)
print("Mean:", MEAN, "Std:", STD)

data = np.nan_to_num(data_nan)  # 增加这一行可以提升DBSCAN速度至少十倍，尚未知晓原因
print("data.shape =", data.shape)

delta_v = np.float64(hdr["CDELT3"]) / 1000
print(delta_v)

cube = sc.read(hdul)
cube = cube.with_spectral_unit(u.km / u.s)
velocities = cube.spectral_axis.value
cube

# data_slab = deepcopy(data)[3478:, :, :]  # 保留[-600,-100]km/s的速度部分

# 合并相邻两个速度切片的值

# stack = data.shape[0]
# for i in range(0, stack, 2):
#     if i // 2 < stack // 2:
#         data[i // 2] = (data[i] + data[i + 1]) / 2
#
# data = data[: stack // 2]

# %%
# hidpi_cube = sc.read('./HI4PI/HI4PI_229_234_4_14.fits')
# hidpi_cube = hidpi_cube.with_spectral_unit(u.km / u.s)
# reprojected_cube = hidpi_cube.reproject(cube.header)
# reprojected_cube.write('./HI4PI/HI4PI_229_234_4_14_reproj.fits', format='fits')

reprojected_cube = sc.read("./HI4PI/HI4PI_229_234_4_14_reproj.fits")
hidpi_data_reproj = reprojected_cube.unmasked_data[:, :, :]

# %%
fig, ax = plt.subplots(1, 4, figsize=[16, 5], subplot_kw={"projection": wcs.celestial})

v = np.arange(100, 225, 25)
for i in range(4):
    slab = cube.spectral_slab(v[i + 1] * u.km / u.s, v[i] * u.km / u.s)
    slab_m0 = slab.moment(order=0)
    norm = simple_norm(np.array(slab_m0), vmin=0, vmax=40, stretch="asinh")
    ax[i].imshow(np.array(slab_m0), norm=norm, cmap="viridis_r", origin="lower")
    ax[i].grid(linestyle="--")
    # lon = ax[i].coords[0]
    # lat = ax[i].coords[1]
    # lon.set_ticks(spacing=5.0 * u.degree)
    # lon.set_axislabel("RA (J2000)")
    # lat.set_axislabel("DEC (J2000)")
    ax[i].set_title(str(v[i]) + "km/s < v < " + str(v[i + 1]) + "km/s")
    ax[i].set_xlabel("Galactic Longitude")
    ax[i].set_ylabel("Galactic Latitude")

fig.suptitle("Moment 0 of Test Area", fontsize=16)
fig.colorbar(mpl.cm.ScalarMappable(cmap=mpl.cm.viridis_r, norm=norm), ax=ax)

# plt.savefig("data.png", dpi=300, bbox_inches='tight')
plt.show()

# %%
# 进行sigma clipping，保留> 3 sigma的点

filtered_data = sigma_clip(data, sigma=3, maxiters=10, masked=True)
mask = filtered_data.mask  # True 的元素为 N > 3 sigma
index = np.transpose(np.nonzero(mask)).astype(int)  # N > 3 sigma的点的坐标[v,y,x]

# %% [markdown]
# ### 记录表
#
# eps = 1, sqrt(2), sqrt(3)
#
# min_samples = 2, 4, 6,
#               6, 9, 12, 15, 18,
#               9, 12, 15, 18
#
# **记录每一个(eps, min_samples)参数组下，表现最好的SNR/spatial_pixels/FWHM**
#
# N为HVC candidate数目，F为其中的假信号数目，T为其中的真信号数目
#
# | eps | minPts | | SNR | PIX | FWHM | | N | F | T |
# | - | - | - | - | - | - | - | - | - | - |
# | 3 | 18 |   | 4 | 6 | 4 |   | 13 | 1 | 12 |
# | 3 | 15 |   | 4 | 6 | 3.5 |   | 14 | 1 | 12 |
# | 3 | 12 |   | 4 | 6 | 3.5 |   | 15 | 2 | 13 |
# | 3 | 9 |    | 3.5 | 6 | 3.5 |   | 14 | 2 | 12 |
# |   |   |    |   |   |   |   |    |   |    |
# | 2 | 18 |   | 4 | 6 | 4 |   | 13 | 1 | 12 |
# | 2 | 15 |   | 4 | 6 | 3.5 |   | 14 | 1 | 13 |
# | 2 | 12 |   | 4 | 6 | 3.5 |   | 15 | 2 | 13 |
# | 2 | 9 |    | 4 | 6 | 4 |   | 12 | 2 | 10 |
# | 2 | 6 |    | 4 | 6 | 4 |   | 12 | 2 | 10 |
# |   |   |    |   |   |   |   |    |   |    |
# | 1 | 6 |    | 4 | 6 | 4 |   | 11 | 1 | 10 |
# | 1 | 4 |    | 4 | 6 | 5 |   | 10 | 1 | 9 |
# | 1 | 2 |    | 4 | 6 | 4 |   | 10 | 1 | 9 |
#
# 备注：eps=3，minPts=9时，(131,70)有一个新的候选体，其余参数下并未见到，可能是被淹没于背景

# %%
# 调整参数

eps = np.sqrt(3)
min_samples = 9

MIN_SNR = 3.5
MIN_PIX = 6
MIN_FWHM = 3.5

# %%
# DBSCAN 算法
# connectivity=1 eps=1 minPts=4
# connectivity=2 eps=sqrt(2) minPts=8
# connectivity=3 eps=sqrt(3) minPts=11

db = DBSCAN(eps=eps, min_samples=min_samples, n_jobs=-1).fit(index)
labels = db.labels_  # 输出的是每一个点属于哪个label

# Number of clusters in labels, ignoring noise if present.
unique_labels = set(labels)
n_clusters = len(unique_labels) - (1 if -1 in labels else 0)
n_noise = list(labels).count(-1)

print("Estimated number of clusters: %d" % n_clusters)
print("Estimated number of noise points: %d" % n_noise)

# %%
# 使用 numpy array 和 spectral cube 计算 moment 0


def moment_0_sc(data, index, labels, n_clusters):
    moment_0_cube = np.empty(0)
    for i in trange(n_clusters):
        vyx = index[labels == i]  # 第i个cluster每个点的坐标
        bool_array = np.zeros(data.shape, dtype=np.bool)
        bool_array[vyx[:, 0], vyx[:, 1], vyx[:, 2]] = (
            True  # 对第i个cluster的每个点标记为True
        )
        bool_cube = sc(bool_array * data, wcs=wcs).with_spectral_unit(u.km / u.s)
        moment_0_coo = sp.coo_array(bool_cube.moment(order=0))
        moment_0_cube = np.hstack((moment_0_cube, moment_0_coo))
    return moment_0_cube


# moment_0_cube_sc = moment_0_sc(data, index, labels, n_clusters)

# %%
# 使用纯 python 循环计算 moment 0


def moment_0_py_vyx(data, vyx):
    # Extract the velocity column and the coordinate column separately
    velo = vyx[:, 0]
    coords = vyx[:, 1:]

    # Find unique rows in the x,y columns and their indices
    unique_coords, indices = np.unique(coords, axis=0, return_inverse=True)
    # print(indices) # 旧列表的元素在新列表的位置

    # Calculate moment 0 of the velocity axis
    # for each unique pair of the coordinates
    moment_0 = np.zeros(len(unique_coords))
    m0 = 0
    for i in range(len(unique_coords)):  # 对每一个坐标
        coo = unique_coords[i]
        for v in velo[indices == i]:  # 对速度轴积分（求和）
            m0 += data[v, coo[0], coo[1]] * abs(delta_v)  # moment 0
        moment_0[i] = m0
        m0 = 0
    moment_0_coo = sp.coo_array(
        (moment_0, (unique_coords[:, 0], unique_coords[:, 1])), shape=data[0].shape
    )
    return moment_0_coo


def moment_0_py(data, index, labels, n_clusters):
    moment_0_cube = np.empty(0)
    for i in trange(n_clusters):
        vyx = index[labels == i]
        moment_0_coo = moment_0_py_vyx(data, vyx)
        moment_0_cube = np.hstack((moment_0_cube, moment_0_coo))
    return moment_0_cube


moment_0_cube = moment_0_py(data, index, labels, n_clusters)

# %%
# 使用 numpy array 和 spectral cube 计算 moment 1


def moment_1_sc(data, index, labels, n_clusters):
    moment_1_cube = np.empty(0)
    for i in trange(n_clusters):
        vyx = index[labels == i]  # 第i个cluster每个点的坐标
        bool_array = np.zeros(data.shape, dtype=np.bool)
        bool_array[vyx[:, 0], vyx[:, 1], vyx[:, 2]] = (
            True  # 对第i个cluster的每个点标记为True
        )
        bool_cube = sc(bool_array * data, wcs=wcs).with_spectral_unit(u.km / u.s)
        moment_1_coo = np.nan_to_num(np.array(bool_cube.moment(order=1)))
        moment_1_coo = sp.coo_array(moment_1_coo)
        moment_1_cube = np.hstack((moment_1_cube, moment_1_coo))
    return np.array(moment_1_cube)


# moment_1_cube_sc = moment_1_sc(data, index, labels, n_clusters)

# %%
# 使用纯 python 循环计算 moment 1


def moment_1_py_vyx(data, vyx):
    # Extract the velocity column and the coordinate column separately
    velo = vyx[:, 0]
    coords = vyx[:, 1:]

    # Find unique rows in the x,y columns and their indices
    unique_coords, indices = np.unique(coords, axis=0, return_inverse=True)
    # print(indices) # 旧列表的元素在新列表的位置

    # Calculate moment 1 of the velocity axis
    # for each unique pair of the coordinates
    moment_1 = np.zeros(len(unique_coords))
    m0 = 0
    m1 = 0
    for i in range(len(unique_coords)):  # 对每一个坐标
        coo = unique_coords[i]
        for v in velo[indices == i]:  # 对速度轴积分（求和）
            m0 += data[v, coo[0], coo[1]]
            vel = velocities[v]
            m1 += data[v, coo[0], coo[1]] * vel  # moment 1
        moment_1[i] = m1 / m0
        m0 = 0
        m1 = 0
    # Combine the unique coordinates with their moment 0
    moment_1_coo = sp.coo_array(
        (moment_1, (unique_coords[:, 0], unique_coords[:, 1])), shape=data[0].shape
    )
    return moment_1_coo


def moment_1_py(data, index, labels, n_clusters):
    moment_1_cube = np.empty(0)
    for i in trange(n_clusters):
        vyx = index[labels == i]
        moment_1_coo = moment_1_py_vyx(data, vyx)
        moment_1_cube = np.hstack((moment_1_cube, moment_1_coo))
    return moment_1_cube


moment_1_cube = moment_1_py(data, index, labels, n_clusters)

# %%
# 使用 numpy array 和 spectral cube 计算 moment 2 (FWHM)


def moment_2_sc(data, index, labels, n_clusters):
    moment_2_cube = np.empty(0)
    for i in trange(n_clusters):
        vyx = index[labels == i]  # 第i个cluster每个点的坐标
        bool_array = np.zeros(data.shape, dtype=np.bool)
        bool_array[vyx[:, 0], vyx[:, 1], vyx[:, 2]] = (
            True  # 对第i个cluster的每个点标记为True
        )
        bool_cube = sc(bool_array * data, wcs=wcs).with_spectral_unit(u.km / u.s)
        moment_2_coo = np.nan_to_num(np.array(bool_cube.linewidth_fwhm()))
        moment_2_coo[np.abs(moment_2_coo) < delta_min] = 0
        moment_2_coo = sp.coo_array(moment_2_coo)
        moment_2_cube = np.hstack((moment_2_cube, moment_2_coo))
    return np.array(moment_2_cube)


# moment_2_cube_sc = moment_2_sc(data, index, labels, n_clusters)

# %%
# 使用纯 python 循环计算 moment 2


def moment_2_py_vyx(data, vyx):
    # Extract the velocity column and the coordinate column separately
    velo = vyx[:, 0]
    coords = vyx[:, 1:]

    # Find unique rows in the x,y columns and their indices
    unique_coords, indices = np.unique(coords, axis=0, return_inverse=True)
    # print(indices) # 旧列表的元素在新列表的位置

    # Calculate moment 1 of the velocity axis
    # for each unique pair of the coordinates
    moment_2 = np.zeros(len(unique_coords))
    m0 = 0
    m1 = 0
    m2 = 0
    for i in range(len(unique_coords)):  # 对每一个坐标
        coo = unique_coords[i]
        for v in velo[indices == i]:  # 对速度轴积分（求和）
            m0 += data[v, coo[0], coo[1]]
            vel = velocities[v]
            m1 += data[v, coo[0], coo[1]] * vel  # moment 1
        m1 = m1 / m0
        for v in velo[indices == i]:
            vel = velocities[v]
            m2 += data[v, coo[0], coo[1]] * (vel - m1) ** 2
        m2 = np.sqrt(8 * np.log(2) * m2 / m0)
        if np.abs(m2) > delta_min:
            moment_2[i] = m2
        else:
            moment_2[i] = 0
        m0 = 0
        m1 = 0
        m2 = 0
    moment_2 = np.nan_to_num(moment_2)
    # Combine the unique coordinates with their moment 0
    moment_2_coo = sp.coo_array(
        (moment_2, (unique_coords[:, 0], unique_coords[:, 1])), shape=data[0].shape
    )
    moment_2_coo.eliminate_zeros()
    return moment_2_coo


def moment_2_py(data, index, labels, n_clusters):
    moment_2_cube = np.empty(0)
    for i in trange(n_clusters):
        vyx = index[labels == i]
        moment_2_coo = moment_2_py_vyx(data, vyx)
        moment_2_cube = np.hstack((moment_2_cube, moment_2_coo))
    return moment_2_cube


moment_2_cube = moment_2_py(data, index, labels, n_clusters)

# %% [markdown]
# # 筛选符合要求的 HVC 候选体
#
# SNR = 3, 4, 5
#
# spatial_pixels = 4, 6, 9  # SIZE = 3.4, 4.2, 5.1 arcmin
#
# FWHM = 3, 4, 5

# %%
# 筛选符合要求的 HVC 候选体

def generate_purity_array(arr, MIN_SNR, MIN_PIX, MIN_FWHM):
    # 使用np.all组合条件
    conditions = np.array(
        [
            arr[:, 0] > MIN_SNR,  # SNR
            arr[:, 1] > MIN_PIX,  # spatial pixels
            (arr[:, 2] > MIN_FWHM) & (arr[:, 2] < 50),  # FWHM
        ]
    )
    TF = np.all(conditions, axis=0)  # 在行方向上检查所有条件
    return TF


def calculate_diameter(spatial_pixels):
    pixel_to_arcmin2 = wcs.proj_plane_pixel_area().to(u.arcmin**2).value
    size = spatial_pixels * pixel_to_arcmin2
    diameter = 2 * np.sqrt(size / np.pi)
    return diameter


def purity(data, index, labels, n_clusters):
    result_array = []
    for i in trange(n_clusters):
        vyx = index[labels == i]  # 第i个cluster每个点的坐标
        clu = data[vyx[:, 0], vyx[:, 1], vyx[:, 2]]
        peak = np.max(clu)
        snr = (peak - MEAN) / STD
        coords = vyx[:, 1:3]
        # 找到唯一的(y, x)组合及其出现次数
        unique_yx, counts = np.unique(coords, axis=0, return_counts=True)
        spatial_pixels = len(unique_yx)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            fwhm = np.nanmean(moment_2_cube[i].data)
        result_array.append([snr, spatial_pixels, fwhm])
    result_array = np.nan_to_num(result_array)
    TF = generate_purity_array(result_array, MIN_SNR, MIN_PIX, MIN_FWHM)
    purity_array = np.tile(TF, (3, 1)).T * result_array
    # 使用np.nonzero函数找到非零行的索引
    real_clusters = np.nonzero(np.any(purity_array != 0, axis=1))[0]
    # 提取非零行
    SNR = purity_array[real_clusters][:, 0]
    SIZE = calculate_diameter(purity_array[real_clusters][:, 1])
    FWHM = purity_array[real_clusters][:, 2]
    return real_clusters, SNR, SIZE, FWHM


real_clusters, SNR, SIZE, FWHM = purity(data, index, labels, n_clusters)
real_clusters

# %%
XCEN = []
YCEN = []
XPEAK = []
YPEAK = []
GLON = []
GLAT = []
RA = []
DEC = []
VLSR = []
TPKB = []
N_HI = []

for i in real_clusters:
    vyx = index[labels == i]  # 第i个cluster每个点的vyx坐标
    bool_array = np.zeros(data.shape, dtype=np.bool_)
    bool_array[vyx[:, 0], vyx[:, 1], vyx[:, 2]] = True
    cluster_data = bool_array * data
    center_vyx = np.mean(vyx, axis=0)
    peak_vyx = np.unravel_index(np.argmax(cluster_data), cluster_data.shape)
    peak_coord = wcs.pixel_to_world_values(peak_vyx[2], peak_vyx[1], peak_vyx[0])
    c_gal = SkyCoord(peak_coord[0], peak_coord[1], frame="galactic", unit="deg")
    c_icrs = c_gal.icrs
    n_HI = np.max(moment_0_cube[i]) * 1.8
    # moment_0_mask = deepcopy(moment_0_cube[i].todense())
    # moment_0_mask[moment_0_mask == 0] = np.nan
    # n_HI = np.nanmean(moment_0_mask) * 1.8
    XCEN.append(center_vyx[2])
    YCEN.append(center_vyx[1])
    XPEAK.append(peak_vyx[2])
    YPEAK.append(peak_vyx[1])
    GLON.append(peak_coord[0])
    GLAT.append(peak_coord[1])
    RA.append(c_icrs.ra.value)
    DEC.append(c_icrs.dec.value)
    VLSR.append(peak_coord[2] / 1000)
    TPKB.append(cluster_data[peak_vyx])
    N_HI.append(n_HI)

XCEN = np.array(XCEN)
YCEN = np.array(YCEN)
XPEAK = np.round(XPEAK) + 1
YPEAK = np.array(YPEAK) + 1
GLON = np.array(GLON)
GLAT = np.array(GLAT)
RA = np.array(RA)
DEC = np.array(DEC)
VLSR = np.array(VLSR)
VGSR = VLSR + 220 * np.sin(GLON) * np.cos(GLAT)
TPKB = np.array(TPKB)
N_HI = np.array(N_HI)

# %%
# SNR, SIZE, FWHM, GLON, GLAT, RA, DEC, VLSR, VGSR, TPKB, N_HI

CATALOG = pd.DataFrame(
    {
        "ID": real_clusters,
        "XPEAK": XPEAK,
        "YPEAK": YPEAK,
        "GLON": GLON,
        "GLAT": GLAT,
        "RA": RA,
        "DEC": DEC,
        "SIZE": SIZE,
        "VLSR": VLSR,
        "VGSR": VGSR,
        "FWHM": FWHM,
        "TPKB": TPKB,
        "SNR": SNR,
    }
)

# CATALOG = CATALOG.sort_values("GLON")
CATALOG

# %%
vyx_clusters = np.empty((0, 3))
for i in real_clusters:
    vyx = index[labels == i]  # 第i个cluster每个点的vyx坐标
    vyx_clusters = np.vstack((vyx_clusters, vyx))
vyx_clusters = vyx_clusters.astype(int)

bool_array = np.zeros(data.shape, dtype=np.bool_)
bool_array[vyx_clusters[:, 0], vyx_clusters[:, 1], vyx_clusters[:, 2]] = (
    True  # 对第i个cluster的每个点标记为True
)
real_clusters_cube = sc(bool_array * data_nan, wcs=wcs).with_spectral_unit(u.km / u.s)

size_pixels = (SIZE / pixel_to_arcmin) ** 2

fig, ax = plt.subplots(
    1, 3, figsize=[12, 6], sharey=True, subplot_kw={"projection": wcs.celestial}
)  # subplot_kw={"projection": wcs.celestial}

moment_0 = np.array(real_clusters_cube.moment(order=0))
im0 = ax[0].imshow(
    moment_0, norm=simple_norm(moment_0, percent=95), cmap="viridis_r", origin="lower"
)
ax[0].scatter(XCEN, YCEN, s=size_pixels, c="None", alpha=1, edgecolor="red")
for i, txt in enumerate(real_clusters):
    ax[0].annotate(i, (XCEN[i], YCEN[i]), va="center", ha="center")
ax[0].grid(linestyle="--")
ax[0].set_title("Moment 0 (K km/s)")
ax[0].set_xlabel("Galactic Longitude")
ax[0].set_ylabel("Galactic Latitude")
plt.colorbar(im0)

moment_1 = np.array(real_clusters_cube.moment(order=1))
im1 = ax[1].imshow(
    moment_1, norm=simple_norm(moment_1, percent=95), cmap="viridis_r", origin="lower"
)
ax[1].scatter(XCEN, YCEN, s=size_pixels, c="None", alpha=1, edgecolor="red")
for i, txt in enumerate(real_clusters):
    ax[1].annotate(i, (XCEN[i], YCEN[i]), va="center", ha="center")
ax[1].grid(linestyle="--")
ax[1].set_title("Moment 1 (km/s)")
ax[1].set_xlabel("Galactic Longitude")
ax[1].set_ylabel("Galactic Latitude")
plt.colorbar(im1)

moment_2 = np.array(real_clusters_cube.linewidth_fwhm())
im2 = ax[2].imshow(
    moment_2, norm=simple_norm(moment_2, percent=95), cmap="viridis_r", origin="lower"
)
ax[2].scatter(XCEN, YCEN, s=size_pixels, c="None", alpha=1, edgecolor="red")
for i, txt in enumerate(real_clusters):
    ax[2].annotate(i, (XCEN[i], YCEN[i]), va="center", ha="center")
ax[2].grid(linestyle="--")
ax[2].set_title("FWHM (km/s)")
ax[2].set_xlabel("Galactic Longitude")
ax[2].set_ylabel("Galactic Latitude")
plt.colorbar(im2)

fig.suptitle(
    "eps = "
    + str(np.round(eps, 3))
    + ", minPts = "
    + str(min_samples)
    + ", MIN_SNR = "
    + str(MIN_SNR)
    + ", MIN_PIX ="
    + str(MIN_PIX)
    + ", FWHM ="
    + str(MIN_FWHM)
)
# plt.savefig("cluster_moments.png", dpi=300, bbox_inches='tight')
plt.show()

# %%
# 对于一个cluster遍历速度轴，对每一个速度v计算所有像素的平均I_v
# 绘制I_v - v光谱图


def calc_spectrum(i, data, index, labels):
    vyx = index[labels == i]  # 第i个cluster每个点的vyx坐标
    bool_array_spec = np.zeros(data.shape, dtype=np.bool_)
    bool_array_spec[:, vyx[:, 1], vyx[:, 2]] = True
    cluster_data = bool_array_spec * data
    cluster_data[cluster_data == 0] = np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        spectrum = np.nanmean(cluster_data, axis=(1, 2))
    vel_array = spectrum != 0
    vel_array = vel_array * velocities
    non_zero_indices = (spectrum != 0) & (vel_array != 0)
    spectrum = spectrum[non_zero_indices]
    vel_array = vel_array[non_zero_indices]
    return vel_array, spectrum

# %%
# HIDPI 的光谱作对比
# 对于一个cluster遍历速度轴，对每一个速度v计算所有像素的平均I_v
# 绘制I_v - v光谱图


vel_array_crafts = []
spectrum_crafts = []

vel_array_hi4pi = []
spectrum_hi4pi = []

for i in real_clusters:
    vel_array1, spectrum1 = calc_spectrum(i, data, index, labels)
    vel_array2, spectrum2 = calc_spectrum(i, hidpi_data_reproj, index, labels)
    vel_array_crafts.append(vel_array1)
    spectrum_crafts.append(spectrum1)
    vel_array_hi4pi.append(vel_array2)
    spectrum_hi4pi.append(spectrum2)

# %% [markdown]
# ### 使用高斯滤波和峰值寻找，根据峰值高度和对称性自动判断光谱信号的真实性

# %%
def get_filter(x_data: np.ndarray, center: float = 0.0, std: float = 1.0, amp: float = 1.0):
    # https://stats.stackexchange.com/a/143633
    return amp * np.exp(-((x_data - center) ** 2) / (2 * std ** 2))

def find_response(x: np.ndarray, y: np.ndarray):
    # Generate the filter you want to match to in the data
    # filter = get_filter(np.linspace(-25, 25, 50), std=3, amp=0.3)
    # y_filtered = np.convolve(y, np.flip(filter), mode='same')
    y_filtered = gaussian_filter(y, sigma=np.mean(FWHM))
    # Use scipy to find all the peaks in the filter response
    peaks, properties = find_peaks(y_filtered)
    # Capture the prominence of each peak (how much does the peak stick up)
    prominences, left_bases, right_bases = peak_prominences(y_filtered, peaks)
    # Find the maximum prominence (which hopefully matches our signal)
    max_prom_index = np.argmax(prominences)
    peak_prominence = prominences[max_prom_index]
    # Find the x-value that the peak occurred at
    main_peak_x = x[peaks[max_prom_index]]
    main_peak_y = y_filtered[peaks[max_prom_index]]
    # the left base of the peak
    main_lbase_x = x[left_bases[max_prom_index]]
    main_lbase_y = y_filtered[left_bases[max_prom_index]]
    # the right base of the peak
    main_rbase_x = x[right_bases[max_prom_index]]
    main_rbase_y = y_filtered[right_bases[max_prom_index]]
    main_peak_plot = [
        main_peak_x,
        main_peak_y,
        main_lbase_x,
        main_lbase_y,
        main_rbase_x,
        main_rbase_y,
    ]
    skew = np.abs(main_peak_plot[2] + main_peak_plot[4] - 2 * main_peak_plot[0])
    skew /= np.abs(main_peak_plot[2] - main_peak_plot[4])
    return y_filtered, main_peak_x, peak_prominence, main_peak_plot, skew

# %%
N_HVC = len(real_clusters)
F_HVC = 0
T_HVC = 0
PEAK_HEIGHT = []
PEAK_WIDTH = []

for i in range(len(vel_array_crafts)):
    y_filtered, v_peak, peak_prominence, main_peak_plot, skew = find_response(
        vel_array_crafts[i], spectrum_crafts[i]
    )
    fig = plt.figure()
    plt.plot(vel_array_crafts[i], spectrum_crafts[i], label="CRAFTS", lw=1)
    plt.plot(vel_array_hi4pi[i], spectrum_hi4pi[i], label="HI4PI", lw=1, ls="--")
    plt.plot(vel_array_crafts[i], y_filtered, label="filtered", lw=1.5)
    # plt.scatter(x[peaks], y_filtered[peaks], c="cyan", marker="+")
    plt.scatter(
        main_peak_plot[0],
        main_peak_plot[1],
        c="cyan",
        marker="D",
        zorder=100,
    )
    # plt.scatter(x[left_bases], y_filtered[left_bases], c="magenta", marker="+")
    plt.scatter(
        main_peak_plot[2],
        main_peak_plot[3],
        c="magenta",
        zorder=100,
    )
    # plt.scatter(x[right_bases], y_filtered[right_bases], c="magenta", marker="+")
    plt.scatter(
        main_peak_plot[4],
        main_peak_plot[5],
        c="magenta",
        zorder=100,
    )
    plt.legend()
    plt.xlim(100, 200)
    plt.ylim(-0.2, 1)
    plt.title("Spectrum of Cluster " + str(i))
    plt.show()
    if peak_prominence > 0.15 and skew < 0.6:
        print("True signal")
        PEAK_HEIGHT.append(peak_prominence)
        PEAK_WIDTH.append(np.abs(main_peak_plot[2] - main_peak_plot[4]))
    else:
        print("False signal")
        F_HVC += 1
    print(i+1, real_clusters[i], v_peak, peak_prominence, skew)


T_HVC = N_HVC - F_HVC
print("N_HVC = " + str(N_HVC) + ", F_HVC = " + str(F_HVC) + ", T_HVC = " + str(T_HVC))
