# %%
import numpy as np

sin = lambda degrees: np.sin(np.deg2rad(degrees))
cos = lambda degrees: np.cos(np.deg2rad(degrees))
tan = lambda degrees: np.tan(np.deg2rad(degrees))

# %%
# r_sun < r_gal < 2 r_sun
def calc_R_min_max_1(l, b, h, r_gal, r_sun):
    l = l % 360
    if b == 0:
        b = 1e-10
    r_b = h / np.abs(tan(b))
    r_l = np.sqrt(r_sun**2 + r_b**2 - 2 * r_sun * r_b * cos(l))
    if r_b >= r_gal + r_sun: # 1
        if cos(l) >= 0:
            R_min = r_sun * np.abs(sin(l))
            R_max = r_gal
        else:
            R_min = r_sun
            R_max = r_gal
    elif 2 * r_sun <= r_b < r_gal + r_sun: # 2
        y_0 = (r_gal**2 + r_sun**2 - r_b**2) / (2 * r_sun)
        if cos(l) >= (r_sun - y_0) / r_b:
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        elif 0 <= cos(l) < (r_sun - y_0) / r_b:
            R_min = r_sun * np.abs(sin(l))
            R_max = r_gal
        else:
            R_min = r_sun
            R_max = r_gal
    elif np.sqrt(r_gal**2 - r_sun**2) <= r_b < 2 * r_sun: # 3
        y_0 = (r_gal**2 + r_sun**2 - r_b**2) / (2 * r_sun)
        if cos(l) >= r_b / (2 * r_sun):
            R_min = r_sun * np.abs(sin(l))
            R_max = r_sun
        elif (r_sun - y_0) / r_b <= cos(l) < r_b / (2 * r_sun):
            R_min = r_sun
            R_max = r_l
        else:
            R_min = r_sun
            R_max = r_gal
    elif r_sun <= r_b < np.sqrt(r_gal**2 - r_sun**2): # 4
        y_0 = (r_gal**2 + r_sun**2 - r_b**2) / (2 * r_sun)
        if cos(l) >= r_b / (2 * r_sun):
            R_min = r_sun * np.abs(sin(l))
            R_max = r_sun
        elif 0 <= cos(l) < r_b / (2 * r_sun):
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        elif (r_sun - y_0) / r_b <= cos(l) < 0:
            R_min = r_sun
            R_max = r_l
        else:
            R_min = r_sun
            R_max = r_gal
    elif r_gal - r_sun <= r_b < r_sun: # 5
        if cos(l) >= r_b / (2 * r_sun):
            R_min = r_l
            R_max = r_sun
        elif r_b / (2 * r_sun) <= cos(l) < r_b / r_sun:
            R_min = r_sun * np.abs(sin(l))
            R_max = r_sun
        elif 0 <= cos(l) < r_b / (2 * r_sun):
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        elif (r_sun - y_0) / r_b <= cos(l) < 0:
            R_min = r_sun
            R_max = r_l
        else:
            R_min = r_sun
            R_max = r_l
    else: # 6
        if cos(l) >= r_b / r_sun:
            R_min = r_l
            R_max = r_sun
        elif r_b / (2 * r_sun) <= cos(l) < r_b / r_sun:
            R_min = r_sun * np.abs(sin(l))
            R_max = r_sun
        elif 0 <= cos(l) < r_b / (2 * r_sun):
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        else:
            R_min = r_sun
            R_max = r_l
    return R_min, R_max

# 2 r_sun < r_gal < sqrt(5) r_sun
def calc_R_min_max_2(l, b, h, r_gal, r_sun):
    l = l % 360
    if b == 0:
        b = 1e-10
    r_b = h / np.abs(tan(b))
    r_l = np.sqrt(r_sun**2 + r_b**2 - 2 * r_sun * r_b * cos(l))
    if r_b >= r_gal + r_sun: # 1
        if cos(l) >= 0:
            R_min = r_sun * np.abs(sin(l))
            R_max = r_gal
        else:
            R_min = r_sun
            R_max = r_gal
    elif 2 * r_sun <= r_b < r_gal + r_sun: # 2
        y_0 = (r_gal**2 + r_sun**2 - r_b**2) / (2 * r_sun)
        if cos(l) >= (r_sun - y_0) / r_b:
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        elif 0 <= cos(l) < (r_sun - y_0) / r_b:
            R_min = r_sun * np.abs(sin(l))
            R_max = r_gal
        else:
            R_min = r_sun
            R_max = r_gal
    elif np.sqrt(r_gal**2 - r_sun**2) <= r_b < 2 * r_sun: # 3
        y_0 = (r_gal**2 + r_sun**2 - r_b**2) / (2 * r_sun)
        if cos(l) >= r_b / (2 * r_sun):
            R_min = r_sun * np.abs(sin(l))
            R_max = r_sun
        elif (r_sun - y_0) / r_b <= cos(l) < r_b / (2 * r_sun):
            R_min = r_sun
            R_max = r_l
        else:
            R_min = r_sun
            R_max = r_gal
    elif r_gal - r_sun <= r_b < np.sqrt(r_gal**2 - r_sun**2): # 4
        y_0 = (r_gal**2 + r_sun**2 - r_b**2) / (2 * r_sun)
        if cos(l) >= r_b / (2 * r_sun):
            R_min = r_sun * np.abs(sin(l))
            R_max = r_sun
        elif 0 <= cos(l) < r_b / (2 * r_sun):
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        elif (r_sun - y_0) / r_b <= cos(l) < 0:
            R_min = r_sun
            R_max = r_l
        else:
            R_min = r_sun
            R_max = r_gal
    elif r_sun <= r_b < r_gal - r_sun: # 5
        if cos(l) >= r_b / (2 * r_sun):
            R_min = r_sun * np.abs(sin(l))
            R_max = r_sun
        elif 0 <= cos(l) < r_b / (2 * r_sun):
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        else:
            R_min = r_sun
            R_max = r_l
    else: # 6
        if cos(l) >= r_b / r_sun:
            R_min = r_l
            R_max = r_sun
        elif r_b / (2 * r_sun) <= cos(l) < r_b / r_sun:
            R_min = r_sun * np.abs(sin(l))
            R_max = r_sun
        elif 0 <= cos(l) < r_b / (2 * r_sun):
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        else:
            R_min = r_sun
            R_max = r_l
    return R_min, R_max

# sqrt(5) r_sun < r_gal < 3 r_sun
def calc_R_min_max_3(l, b, h, r_gal, r_sun):
    l = l % 360
    if b == 0:
        b = 1e-10
    r_b = h / np.abs(tan(b))
    r_l = np.sqrt(r_sun**2 + r_b**2 - 2 * r_sun * r_b * cos(l))
    if r_b >= r_gal + r_sun: # 1
        if cos(l) >= 0:
            R_min = r_sun * np.abs(sin(l))
            R_max = r_gal
        else:
            R_min = r_sun
            R_max = r_gal
    elif np.sqrt(r_gal**2 - r_sun**2) <= r_b < r_gal + r_sun: # 2
        y_0 = (r_gal**2 + r_sun**2 - r_b**2) / (2 * r_sun)
        if cos(l) >= (r_sun - y_0) / r_b:
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        elif 0 <= cos(l) < (r_sun - y_0) / r_b:
            R_min = r_sun * np.abs(sin(l))
            R_max = r_gal
        else:
            R_min = r_sun
            R_max = r_gal
    elif 2 * r_sun <= r_b < np.sqrt(r_gal**2 - r_sun**2): # 3
        y_0 = (r_gal**2 + r_sun**2 - r_b**2) / (2 * r_sun)
        if cos(l) >= 0:
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        elif (r_sun - y_0) / r_b <= cos(l) < 0:
            R_min = r_sun
            R_max = r_l
        else:
            R_min = r_sun
            R_max = r_gal
    elif r_gal - r_sun <= r_b < 2 * r_sun: # 4
        y_0 = (r_gal**2 + r_sun**2 - r_b**2) / (2 * r_sun)
        if cos(l) >= r_b / (2 * r_sun):
            R_min = r_sun * np.abs(sin(l))
            R_max = r_sun
        elif 0 <= cos(l) < r_b / (2 * r_sun):
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        elif (r_sun - y_0) / r_b <= cos(l) < 0:
            R_min = r_sun
            R_max = r_l
        else:
            R_min = r_sun
            R_max = r_gal
    elif r_sun <= r_b < r_gal - r_sun: # 5
        if cos(l) >= r_b / (2 * r_sun):
            R_min = r_sun * np.abs(sin(l))
            R_max = r_sun
        elif 0 <= cos(l) < r_b / (2 * r_sun):
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        else:
            R_min = r_sun
            R_max = r_l
    else: # 6
        if cos(l) >= r_b / r_sun:
            R_min = r_l
            R_max = r_sun
        elif r_b / (2 * r_sun) <= cos(l) < r_b / r_sun:
            R_min = r_sun * np.abs(sin(l))
            R_max = r_sun
        elif 0 <= cos(l) < r_b / (2 * r_sun):
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        else:
            R_min = r_sun
            R_max = r_l
    return R_min, R_max

# r_gal >= 3 r_sun
def calc_R_min_max_4(l, b, h, r_gal, r_sun):
    l = l % 360
    if b == 0:
        b = 1e-10
    r_b = h / np.abs(tan(b))
    r_l = np.sqrt(r_sun**2 + r_b**2 - 2 * r_sun * r_b * cos(l))
    if r_b >= r_gal + r_sun: # 1
        if cos(l) >= 0:
            R_min = r_sun * np.abs(sin(l))
            R_max = r_gal
        else:
            R_min = r_sun
            R_max = r_gal
    elif np.sqrt(r_gal**2 - r_sun**2) <= r_b < r_gal + r_sun: # 2
        y_0 = (r_gal**2 + r_sun**2 - r_b**2) / (2 * r_sun)
        if cos(l) >= (r_sun - y_0) / r_b:
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        elif 0 <= cos(l) < (r_sun - y_0) / r_b:
            R_min = r_sun * np.abs(sin(l))
            R_max = r_gal
        else:
            R_min = r_sun
            R_max = r_gal
    elif r_gal - r_sun <= r_b < np.sqrt(r_gal**2 - r_sun**2): # 3
        y_0 = (r_gal**2 + r_sun**2 - r_b**2) / (2 * r_sun)
        if cos(l) >= 0:
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        elif (r_sun - y_0) / r_b <= cos(l) < 0:
            R_min = r_sun
            R_max = r_l
        else:
            R_min = r_sun
            R_max = r_gal
    elif 2 * r_sun <= r_b < r_gal - r_sun: # 4
        if cos(l) >= 0:
            R_min = r_sun * np.abs(sin(l))
            R_max = r_sun
        else:
            R_min = r_sun
            R_max = r_l
    elif r_sun <= r_b < 2 * r_sun: # 5
        if cos(l) >= r_b / (2 * r_sun):
            R_min = r_sun * np.abs(sin(l))
            R_max = r_sun
        elif 0 <= cos(l) < r_b / (2 * r_sun):
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        else:
            R_min = r_sun
            R_max = r_l
    else: # 6
        if cos(l) >= r_b / r_sun:
            R_min = r_l
            R_max = r_sun
        elif r_b / (2 * r_sun) <= cos(l) < r_b / r_sun:
            R_min = r_sun * np.abs(sin(l))
            R_max = r_sun
        elif 0 <= cos(l) < r_b / (2 * r_sun):
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        else:
            R_min = r_sun
            R_max = r_l
    return R_min, R_max


# %%
# Simple rotation curve model
# V_rot = 220 km/s where R > 0.5 kpc, else V_rot = 440 * R
# Wakker 1991: 1991A&A...250..499W
def v_rot_simple(r, v_sun=220, r_cut=0.5):
    v_rot = np.where(r > r_cut, v_sun, (1 / r_cut) * v_sun * r)
    return v_rot

# %%
# Universal rotation curve model
# Persic et al. 1996: 1996MNRAS.281...27P
# Reid et al. 2019: 2019ApJ...885..131R
def v_rot_univ(r, a2=0.96, a3=1.62):
    """
    Calculate the circular rotation speed (Vr) at a given galactocentric radius (r).

    Parameters:
        r (float): Galactocentric radius in kpc.
        a2 (float): r_opt/r_sun, where r_opt = 3.2 * r_scale_length encloses 83% of light
        a3 (float): 1.5 * (L/L*)^(1/5)
        r_sun (float): The distance from Sun to Galactic Center in kpc.

    Returns:
        Vr (float): Circular rotation speed at radius r in km/s.
    """
    # Calculate lambda = L/L*
    lambda_ = (a3 / 1.5) ** 5
    # Calculate r_opt
    r_sun = 8.15
    r_opt = a2 * r_sun
    # Calculate normalized radius x
    x = r / r_opt
    # Calculate log10(lambda)
    log_lam = np.log10(lambda_)
    # Calculate term1
    term1 = 200.0 * lambda_**0.41
    # Calculate term2
    top = 0.75 * np.exp(-0.4 * lambda_)
    bot = 0.47 + 2.25 * lambda_**0.4
    term2 = np.sqrt(0.80 + 0.49 * log_lam + (top / bot))
    # Calculate term3
    top = 1.97 * x**1.22
    bot = (x**2 + 0.61) ** 1.43
    term3 = (0.72 + 0.44 * log_lam) * (top / bot)
    # Calculate term4
    top = x**2
    bot = x**2 + 2.25 * lambda_**0.4
    term4 = 1.6 * np.exp(-0.4 * lambda_) * (top / bot)
    # Calculate Vr (rotation speed)
    Vr = (term1 / term2) * np.sqrt(term3 + term4)  # km/s
    return Vr

# %%
# Linear rotation curve model
# V_rot = 229 - 1.7 * (R - R_sun)
# Eilers et al. 2019: 2019ApJ...871...120E
def v_rot_linear(r):
    r_sun = 8.122
    v_rot = 229 - 1.7 * (r - r_sun)
    return v_rot

# %%
# Power law rotation curve model
# V_rot = 240 * 1.022 * (R / 8.34)^0.0803
# Russeil et al. 2017: 2017A&A...601L...5R
def v_rot_power(r):
    r_sun = 8.34
    v_sun = 240
    v_rot = v_sun * 1.022 * (r / r_sun)**0.0803
    return v_rot

# %%
# Calculate V_max and V_min for given (l, b, h, r_gal)
# Simple rotation curve model
def calc_v_max_min_simple(l, b, h, r_gal, r_sun=8.5, v_sun = 220):
    v_sun = v_rot_simple(r_sun, v_sun)
    if r_sun < r_gal < 2 * r_sun:
        R_min, R_max = calc_R_min_max_1(l, b, h, r_gal, r_sun)
    elif 2 * r_sun < r_gal < np.sqrt(5) * r_sun:
        R_min, R_max = calc_R_min_max_2(l, b, h, r_gal, r_sun)
    elif np.sqrt(5) * r_sun < r_gal < 3 * r_sun:
        R_min, R_max = calc_R_min_max_3(l, b, h, r_gal, r_sun)
    else:
        R_min, R_max = calc_R_min_max_4(l, b, h, r_gal, r_sun)
    v_R_min = (v_rot_simple(R_min) * r_sun / R_min - v_sun) * sin(l) * cos(b)
    v_R_max = (v_rot_simple(R_max) * r_sun / R_max - v_sun) * sin(l) * cos(b)
    v_max = np.maximum(v_R_min, v_R_max)
    v_min = np.minimum(v_R_min, v_R_max)
    return v_max, v_min


# %%
# Calculate V_max and V_min for given (l, b, h, r_gal)
# Universal rotation curve model
def calc_v_max_min_univ(l, b, h, r_gal):
    r_sun = 8.15
    v_sun = v_rot_univ(r_sun)
    if r_sun < r_gal < 2 * r_sun:
        R_min, R_max = calc_R_min_max_1(l, b, h, r_gal, r_sun)
    elif 2 * r_sun < r_gal < np.sqrt(5) * r_sun:
        R_min, R_max = calc_R_min_max_2(l, b, h, r_gal, r_sun)
    elif np.sqrt(5) * r_sun < r_gal < 3 * r_sun:
        R_min, R_max = calc_R_min_max_3(l, b, h, r_gal, r_sun)
    else:
        R_min, R_max = calc_R_min_max_4(l, b, h, r_gal, r_sun)
    v_R_min = (v_rot_univ(R_min) * r_sun / R_min - v_sun) * sin(l) * cos(b)
    v_R_max = (v_rot_univ(R_max) * r_sun / R_max - v_sun) * sin(l) * cos(b)
    v_max = np.maximum(v_R_min, v_R_max)
    v_min = np.minimum(v_R_min, v_R_max)
    return v_max, v_min


# %%
# Calculate V_max and V_min for given (l, b, h, r_gal)
# Linear rotation curve model
def calc_v_max_min_linear(l, b, h, r_gal):
    r_sun = 8.122
    v_sun = v_rot_linear(r_sun)
    if r_sun < r_gal < 2 * r_sun:
        R_min, R_max = calc_R_min_max_1(l, b, h, r_gal, r_sun)
    elif 2 * r_sun < r_gal < np.sqrt(5) * r_sun:
        R_min, R_max = calc_R_min_max_2(l, b, h, r_gal, r_sun)
    elif np.sqrt(5) * r_sun < r_gal < 3 * r_sun:
        R_min, R_max = calc_R_min_max_3(l, b, h, r_gal, r_sun)
    else:
        R_min, R_max = calc_R_min_max_4(l, b, h, r_gal, r_sun)
    if R_min == 0:
        R_min = 1e-10
    if R_max == 0:
        R_max = 1e-10
    v_R_min = (v_rot_linear(R_min) * r_sun / R_min - v_sun) * sin(l) * cos(b)
    v_R_max = (v_rot_linear(R_max) * r_sun / R_max - v_sun) * sin(l) * cos(b)
    v_max = np.maximum(v_R_min, v_R_max)
    v_min = np.minimum(v_R_min, v_R_max)
    return v_max, v_min


# %%
# Calculate V_max and V_min for given (l, b, h, r_gal)
# Power law rotation curve model
def calc_v_max_min_power(l, b, h, r_gal, v_sun = 220):
    r_sun = 8.34
    v_sun = v_rot_power(r_sun)
    if r_sun < r_gal < 2 * r_sun:
        R_min, R_max = calc_R_min_max_1(l, b, h, r_gal, r_sun)
    elif 2 * r_sun < r_gal < np.sqrt(5) * r_sun:
        R_min, R_max = calc_R_min_max_2(l, b, h, r_gal, r_sun)
    elif np.sqrt(5) * r_sun < r_gal < 3 * r_sun:
        R_min, R_max = calc_R_min_max_3(l, b, h, r_gal, r_sun)
    else:
        R_min, R_max = calc_R_min_max_4(l, b, h, r_gal, r_sun)
    v_R_min = (v_rot_power(R_min) * r_sun / R_min - v_sun) * sin(l) * cos(b)
    v_R_max = (v_rot_power(R_max) * r_sun / R_max - v_sun) * sin(l) * cos(b)
    v_max = np.maximum(v_R_min, v_R_max)
    v_min = np.minimum(v_R_min, v_R_max)
    return v_max, v_min


# %%
# Deviation Velocity
def calc_v_dev_py(l, b, h=5, r_gal=20, r_sun=8.5, v_sun=220, model='univ', v_dev=50):
    if model == 'simple':
        v_max, v_min = calc_v_max_min_simple(l, b, h, r_gal, r_sun, v_sun)
    elif model == 'univ':
        v_max, v_min = calc_v_max_min_univ(l, b, h, r_gal)
    elif model == 'linear':
        v_max, v_min = calc_v_max_min_linear(l, b, h, r_gal)
    elif model == 'power':
        v_max, v_min = calc_v_max_min_power(l, b, h, r_gal)
    else:
        raise ValueError("Invalid model. Choose from 'simple', 'univ', 'linear', 'power'.")
    return v_max + v_dev, v_min - v_dev