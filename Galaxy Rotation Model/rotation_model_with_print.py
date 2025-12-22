import numpy as np
from numpy.polynomial import Polynomial as P

# Degree to Radian conversion trigonometric functions
sin = lambda degrees: np.sin(np.deg2rad(degrees))
cos = lambda degrees: np.cos(np.deg2rad(degrees))
tan = lambda degrees: np.tan(np.deg2rad(degrees))

# r_sun < r_gal < 2 r_sun
def calc_R_min_max_1(l, b, h, r_gal, r_sun):
    print("1. r_sun < r_gal < 2 r_sun")
    l = l % 360
    if b == 0:
        b = 1e-10
    r_b = h / np.abs(tan(b))
    r_l = np.sqrt(r_sun**2 + r_b**2 - 2 * r_sun * r_b * cos(l))
    y_0 = (r_gal**2 + r_sun**2 - r_b**2) / (2 * r_sun)
    if r_b >= r_gal + r_sun:  # 1
        print("condition 1: r_b >= r_gal + r_sun")
        if cos(l) >= 0:
            print("sub-condition 1.1: cos(l) >= 0")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_gal
        else:
            print("sub-condition 1.2: cos(l) < 0")
            R_min = r_sun
            R_max = r_gal
    elif 2 * r_sun <= r_b < r_gal + r_sun:  # 2
        print("condition 2: 2 * r_sun <= r_b < r_gal + r_sun")
        if cos(l) >= (r_sun - y_0) / r_b:
            print("sub-condition 2.1: cos(l) >= (r_sun - y_0) / r_b")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        elif 0 <= cos(l) < (r_sun - y_0) / r_b:
            print("sub-condition 2.2: 0 <= cos(l) < (r_sun - y_0) / r_b")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_gal
        else:
            print("sub-condition 2.3: cos(l) < 0")
            R_min = r_sun
            R_max = r_gal
    elif np.sqrt(r_gal**2 - r_sun**2) <= r_b < 2 * r_sun:  # 3
        print("condition 3: np.sqrt(r_gal**2 - r_sun**2) <= r_b < 2 * r_sun")
        if cos(l) >= r_b / (2 * r_sun):
            print("sub-condition 3.1: cos(l) >= r_b / (2 * r_sun)")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_sun
        elif (r_sun - y_0) / r_b <= cos(l) < r_b / (2 * r_sun):
            print("sub-condition 3.2: (r_sun - y_0) / r_b <= cos(l) < r_b / (2 * r_sun)")
            R_min = r_sun
            R_max = r_l
        else:
            print("sub-condition 3.3: cos(l) < (r_sun - y_0) / r_b")
            R_min = r_sun
            R_max = r_gal
    elif r_sun <= r_b < np.sqrt(r_gal**2 - r_sun**2):  # 4
        print("condition 4: r_sun <= r_b < np.sqrt(r_gal**2 - r_sun**2)")
        if cos(l) >= r_b / (2 * r_sun):
            print("sub-condition 4.1: cos(l) >= r_b / (2 * r_sun)")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_sun
        elif 0 <= cos(l) < r_b / (2 * r_sun):
            print("sub-condition 4.2: 0 <= cos(l) < r_b / (2 * r_sun)")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        elif (r_sun - y_0) / r_b <= cos(l) < 0:
            print("sub-condition 4.3: (r_sun - y_0) / r_b <= cos(l) < 0")
            R_min = r_sun
            R_max = r_l
        else:
            print("sub-condition 4.4: cos(l) < (r_sun - y_0) / r_b")
            R_min = r_sun
            R_max = r_gal
    elif r_gal - r_sun <= r_b < r_sun:  # 5
        print("condition 5: r_gal - r_sun <= r_b < r_sun")
        if cos(l) >= r_b / r_sun:
            print("sub-condition 5.1: cos(l) >= r_b / r_sun")
            R_min = r_l
            R_max = r_sun
        elif r_b / (2 * r_sun) <= cos(l) < r_b / r_sun:
            print("sub-condition 5.2: r_b / (2 * r_sun) <= cos(l) < r_b / r_sun")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_sun
        elif 0 <= cos(l) < r_b / (2 * r_sun):
            print("sub-condition 5.3: 0 <= cos(l) < r_b / (2 * r_sun)")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        elif (r_sun - y_0) / r_b <= cos(l) < 0:
            print("sub-condition 5.4: (r_sun - y_0) / r_b <= cos(l) < 0")
            R_min = r_sun
            R_max = r_l
        else:
            print("sub-condition 5.5: cos(l) < (r_sun - y_0) / r_b")
            R_min = r_sun
            R_max = r_l
    else:  # 6
        print("condition 6: r_b < r_gal - r_sun")
        if cos(l) >= r_b / r_sun:
            print("sub-condition 6.1: cos(l) >= r_b / r_sun")
            R_min = r_l
            R_max = r_sun
        elif r_b / (2 * r_sun) <= cos(l) < r_b / r_sun:
            print("sub-condition 6.2: r_b / (2 * r_sun) <= cos(l) < r_b / r_sun")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_sun
        elif 0 <= cos(l) < r_b / (2 * r_sun):
            print("sub-condition 6.3: 0 <= cos(l) < r_b / (2 * r_sun)")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        else:
            print("sub-condition 6.4: cos(l) < 0")
            R_min = r_sun
            R_max = r_l
    return R_min, R_max

# 2 r_sun < r_gal < sqrt(5) r_sun
def calc_R_min_max_2(l, b, h, r_gal, r_sun):
    print("2. 2 r_sun < r_gal < sqrt(5) r_sun")
    l = l % 360
    if b == 0:
        b = 1e-10
    r_b = h / np.abs(tan(b))
    r_l = np.sqrt(r_sun**2 + r_b**2 - 2 * r_sun * r_b * cos(l))
    y_0 = (r_gal**2 + r_sun**2 - r_b**2) / (2 * r_sun)
    if r_b >= r_gal + r_sun:  # 1
        print("condition 1: r_b >= r_gal + r_sun")
        if cos(l) >= 0:
            print("sub-condition 1.1: cos(l) >= 0")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_gal
        else:
            print("sub-condition 1.2: cos(l) < 0")
            R_min = r_sun
            R_max = r_gal
    elif 2 * r_sun <= r_b < r_gal + r_sun:  # 2
        print("condition 2: 2 * r_sun <= r_b < r_gal + r_sun")
        if cos(l) >= (r_sun - y_0) / r_b:
            print("sub-condition 2.1: cos(l) >= (r_sun - y_0) / r_b")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        elif 0 <= cos(l) < (r_sun - y_0) / r_b:
            print("sub-condition 2.2: 0 <= cos(l) < (r_sun - y_0) / r_b")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_gal
        else:
            print("sub-condition 2.3: cos(l) < 0")
            R_min = r_sun
            R_max = r_gal
    elif np.sqrt(r_gal**2 - r_sun**2) <= r_b < 2 * r_sun:  # 3
        print("condition 3: np.sqrt(r_gal**2 - r_sun**2) <= r_b < 2 * r_sun")
        if cos(l) >= r_b / (2 * r_sun):
            print("sub-condition 3.1: cos(l) >= r_b / (2 * r_sun)")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_sun
        elif (r_sun - y_0) / r_b <= cos(l) < r_b / (2 * r_sun):
            print("sub-condition 3.2: (r_sun - y_0) / r_b <= cos(l) < r_b / (2 * r_sun)")
            R_min = r_sun
            R_max = r_l
        else:
            print("sub-condition 3.3: cos(l) < (r_sun - y_0) / r_b")
            R_min = r_sun
            R_max = r_gal
    elif r_gal - r_sun <= r_b < np.sqrt(r_gal**2 - r_sun**2):  # 4
        print("condition 4: r_gal - r_sun <= r_b < np.sqrt(r_gal**2 - r_sun**2)")
        if cos(l) >= r_b / (2 * r_sun):
            print("sub-condition 4.1: cos(l) >= r_b / (2 * r_sun)")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_sun
        elif 0 <= cos(l) < r_b / (2 * r_sun):
            print("sub-condition 4.2: 0 <= cos(l) < r_b / (2 * r_sun)")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        elif (r_sun - y_0) / r_b <= cos(l) < 0:
            print("sub-condition 4.3: (r_sun - y_0) / r_b <= cos(l) < 0")
            R_min = r_sun
            R_max = r_l
        else:
            print("sub-condition 4.4: cos(l) < (r_sun - y_0) / r_b")
            R_min = r_sun
            R_max = r_gal
    elif r_sun <= r_b < r_gal - r_sun:  # 5
        print("condition 5: r_sun <= r_b < r_gal - r_sun")
        if cos(l) >= r_b / (2 * r_sun):
            print("sub-condition 5.1: cos(l) >= r_b / (2 * r_sun)")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_sun
        elif 0 <= cos(l) < r_b / (2 * r_sun):
            print("sub-condition 5.2: 0 <= cos(l) < r_b / (2 * r_sun)")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        else:
            print("sub-condition 5.3: cos(l) < 0")
            R_min = r_sun
            R_max = r_l
    else:  # 6
        print("condition 6: r_b < r_sun")
        if cos(l) >= r_b / r_sun:
            print("sub-condition 6.1: cos(l) >= r_b / r_sun")
            R_min = r_l
            R_max = r_sun
        elif r_b / (2 * r_sun) <= cos(l) < r_b / r_sun:
            print("sub-condition 6.2: r_b / (2 * r_sun) <= cos(l) < r_b / r_sun")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_sun
        elif 0 <= cos(l) < r_b / (2 * r_sun):
            print("sub-condition 6.3: 0 <= cos(l) < r_b / (2 * r_sun)")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        else:
            print("sub-condition 6.4: cos(l) < 0")
            R_min = r_sun
            R_max = r_l
    return R_min, R_max

# sqrt(5) r_sun < r_gal < 3 r_sun
def calc_R_min_max_3(l, b, h, r_gal, r_sun):
    print("3. sqrt(5) r_sun < r_gal < 3 r_sun")
    l = l % 360
    if b == 0:
        b = 1e-10
    r_b = h / np.abs(tan(b))
    r_l = np.sqrt(r_sun**2 + r_b**2 - 2 * r_sun * r_b * cos(l))
    y_0 = (r_gal**2 + r_sun**2 - r_b**2) / (2 * r_sun)
    if r_b >= r_gal + r_sun:  # 1
        print("condition 1: r_b >= r_gal + r_sun")
        if cos(l) >= 0:
            print("sub-condition 1.1: cos(l) >= 0")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_gal
        else:
            print("sub-condition 1.2: cos(l) < 0")
            R_min = r_sun
            R_max = r_gal
    elif np.sqrt(r_gal**2 - r_sun**2) <= r_b < r_gal + r_sun:  # 2
        print("condition 2: np.sqrt(r_gal**2 - r_sun**2) <= r_b < r_gal + r_sun")
        if cos(l) >= (r_sun - y_0) / r_b:
            print("sub-condition 2.1: cos(l) >= (r_sun - y_0) / r_b")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        elif 0 <= cos(l) < (r_sun - y_0) / r_b:
            print("sub-condition 2.2: 0 <= cos(l) < (r_sun - y_0) / r_b")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_gal
        else:
            print("sub-condition 2.3: cos(l) < 0")
            R_min = r_sun
            R_max = r_gal
    elif 2 * r_sun <= r_b < np.sqrt(r_gal**2 - r_sun**2):  # 3
        print("condition 3: 2 * r_sun <= r_b < np.sqrt(r_gal**2 - r_sun**2)")
        if cos(l) >= 0:
            print("sub-condition 3.1: cos(l) >= 0")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        elif (r_sun - y_0) / r_b <= cos(l) < 0:
            print("sub-condition 3.2: (r_sun - y_0) / r_b <= cos(l) < 0")
            R_min = r_sun
            R_max = r_l
        else:
            print("sub-condition 3.3: cos(l) < (r_sun - y_0) / r_b")
            R_min = r_sun
            R_max = r_gal
    elif r_gal - r_sun <= r_b < 2 * r_sun:  # 4
        print("condition 4: r_gal - r_sun <= r_b < 2 * r_sun")
        if cos(l) >= r_b / (2 * r_sun):
            print("sub-condition 4.1: cos(l) >= r_b / (2 * r_sun)")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_sun
        elif 0 <= cos(l) < r_b / (2 * r_sun):
            print("sub-condition 4.2: 0 <= cos(l) < r_b / (2 * r_sun)")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        elif (r_sun - y_0) / r_b <= cos(l) < 0:
            print("sub-condition 4.3: (r_sun - y_0) / r_b <= cos(l) < 0")
            R_min = r_sun
            R_max = r_l
        else:
            print("sub-condition 4.4: cos(l) < (r_sun - y_0) / r_b")
            R_min = r_sun
            R_max = r_gal
    elif r_sun <= r_b < r_gal - r_sun:  # 5
        print("condition 5: r_sun <= r_b < r_gal - r_sun")
        if cos(l) >= r_b / (2 * r_sun):
            print("sub-condition 5.1: cos(l) >= r_b / (2 * r_sun)")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_sun
        elif 0 <= cos(l) < r_b / (2 * r_sun):
            print("sub-condition 5.2: 0 <= cos(l) < r_b / (2 * r_sun)")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        else:
            print("sub-condition 5.3: cos(l) < 0")
            R_min = r_sun
            R_max = r_l
    else:  # 6
        print("condition 6: r_b < r_sun")
        if cos(l) >= r_b / r_sun:
            print("sub-condition 6.1: cos(l) >= r_b / r_sun")
            R_min = r_l
            R_max = r_sun
        elif r_b / (2 * r_sun) <= cos(l) < r_b / r_sun:
            print("sub-condition 6.2: r_b / (2 * r_sun) <= cos(l) < r_b / r_sun")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_sun
        elif 0 <= cos(l) < r_b / (2 * r_sun):
            print("sub-condition 6.3: 0 <= cos(l) < r_b / (2 * r_sun)")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        else:
            print("sub-condition 6.4: cos(l) < 0")
            R_min = r_sun
            R_max = r_l
    return R_min, R_max

# r_gal >= 3 r_sun
def calc_R_min_max_4(l, b, h, r_gal, r_sun):
    print("4. r_gal >= 3 r_sun")
    l = l % 360
    if b == 0:
        b = 1e-10
    r_b = h / np.abs(tan(b))
    r_l = np.sqrt(r_sun**2 + r_b**2 - 2 * r_sun * r_b * cos(l))
    y_0 = (r_gal**2 + r_sun**2 - r_b**2) / (2 * r_sun)
    if r_b >= r_gal + r_sun:  # 1
        print("condition 1: r_b >= r_gal + r_sun")
        if cos(l) >= 0:
            print("sub-condition 1.1: cos(l) >= 0")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_gal
        else:
            print("sub-condition 1.2: cos(l) < 0")
            R_min = r_sun
            R_max = r_gal
    elif np.sqrt(r_gal**2 - r_sun**2) <= r_b < r_gal + r_sun:  # 2
        print("condition 2: np.sqrt(r_gal**2 - r_sun**2) <= r_b < r_gal + r_sun")
        if cos(l) >= (r_sun - y_0) / r_b:
            print("sub-condition 2.1: cos(l) >= (r_sun - y_0) / r_b")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        elif 0 <= cos(l) < (r_sun - y_0) / r_b:
            print("sub-condition 2.2: 0 <= cos(l) < (r_sun - y_0) / r_b")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_gal
        else:
            print("sub-condition 2.3: cos(l) < 0")
            R_min = r_sun
            R_max = r_gal
    elif r_gal - r_sun <= r_b < np.sqrt(r_gal**2 - r_sun**2):  # 3
        print("condition 3: r_gal - r_sun <= r_b < np.sqrt(r_gal**2 - r_sun**2)")
        if cos(l) >= 0:
            print("sub-condition 3.1: cos(l) >= 0")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        elif (r_sun - y_0) / r_b <= cos(l) < 0:
            print("sub-condition 3.2: (r_sun - y_0) / r_b <= cos(l) < 0")
            R_min = r_sun
            R_max = r_l
        else:
            print("sub-condition 3.3: cos(l) < (r_sun - y_0) / r_b")
            R_min = r_sun
            R_max = r_gal
    elif 2 * r_sun <= r_b < r_gal - r_sun:  # 4
        print("condition 4: 2 * r_sun <= r_b < r_gal - r_sun")
        if cos(l) >= 0:
            print("sub-condition 4.1: cos(l) >= 0")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_sun
        else:
            print("sub-condition 4.2: cos(l) < 0")
            R_min = r_sun
            R_max = r_l
    elif r_sun <= r_b < 2 * r_sun:  # 5
        print("condition 5: r_sun <= r_b < 2 * r_sun")
        if cos(l) >= r_b / (2 * r_sun):
            print("sub-condition 5.1: cos(l) >= r_b / (2 * r_sun)")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_sun
        elif 0 <= cos(l) < r_b / (2 * r_sun):
            print("sub-condition 5.2: 0 <= cos(l) < r_b / (2 * r_sun)")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        else:
            print("sub-condition 5.3: cos(l) < 0")
            R_min = r_sun
            R_max = r_l
    else:  # 6
        print("condition 6: r_b < r_sun")
        if cos(l) >= r_b / r_sun:
            print("sub-condition 6.1: cos(l) >= r_b / r_sun")
            R_min = r_l
            R_max = r_sun
        elif r_b / (2 * r_sun) <= cos(l) < r_b / r_sun:
            print("sub-condition 6.2: r_b / (2 * r_sun) <= cos(l) < r_b / r_sun")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_sun
        elif 0 <= cos(l) < r_b / (2 * r_sun):
            print("sub-condition 6.3: 0 <= cos(l) < r_b / (2 * r_sun)")
            R_min = r_sun * np.abs(sin(l))
            R_max = r_l
        else:
            print("sub-condition 6.4: cos(l) < 0")
            R_min = r_sun
            R_max = r_l
    return R_min, R_max

# Simple rotation curve model
# V_rot = 220 km/s where R > 0.5 kpc, else V_rot = 440 * R
# Wakker 1991: 1991A&A...250..499W
def v_rot_simple(r, v_sun=220, r_cut=0.5):
    v_rot = np.where(r > r_cut, v_sun, (1 / r_cut) * v_sun * r)
    return v_rot

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

# Linear rotation curve model
# V_rot = 234.04 - 1.83 * (R - R_sun)
# Zhou et al. 2023: 2023ApJ...946...73Z
def v_rot_linear(r):
    r_sun = 8.122
    v_rot = 234.04 - 1.83 * (r - r_sun)
    return v_rot

# Combined Polynomial rotation curve model
# Clemens 1985: 1985ApJ...295..422C
def v_rot_poly(r):
    v_rot = np.zeros_like(r)
    r_sun = 8.5
    # conditions for r
    condition1 = r < 0.09 * r_sun
    condition2 = (0.09 * r_sun <= r) & (r < 0.45 * r_sun)
    condition3 = (0.45 * r_sun <= r) & (r < 1.6 * r_sun)
    condition4 = r >= 1.6 * r_sun
    # r < 0.09 * r_sun
    coeffs_A = [0, 3069.81, -15809.8, 43980.1, -68287.3, 54904, -17731]
    poly_A = P(coeffs_A)
    v_rot = np.where(condition1, poly_A(r), v_rot)
    # (0.09 * r_sun <= r) & (r < 0.45 * r_sun)
    coeffs_B = [325.0912, -248.1467, 231.87099, -110.73531, 25.073006, -2.110625]
    poly_B = P(coeffs_B)
    v_rot = np.where(condition2, poly_B(r), v_rot)
    # (0.45 * r_sun <= r) & (r < 1.6 * r_sun)
    coeffs_C = [
        -2342.6564,
        2507.60391,
        -1024.068760,
        224.562732,
        -28.4080026,
        2.0697271,
        -0.08050808,
        0.00129348,
    ]
    poly_C = P(coeffs_C)
    v_rot = np.where(condition3, poly_C(r), v_rot)
    # r >= 1.6 * r_sun
    v_rot = np.where(condition4, 234.88, v_rot)
    return v_rot

# Calculate V_max and V_min for given (l, b, h, r_gal)
# Simple rotation curve model
def calc_v_max_min_simple(l, b, h, r_gal, r_sun=8.5, v_sun=220, r_cut=0.5):
    v_sun = v_rot_simple(r_sun, v_sun, r_cut)
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
    v_R_min = (v_rot_simple(R_min) * r_sun / R_min - v_sun) * sin(l) * cos(b)
    v_R_max = (v_rot_simple(R_max) * r_sun / R_max - v_sun) * sin(l) * cos(b)
    v_max = np.maximum(v_R_min, v_R_max)
    v_min = np.minimum(v_R_min, v_R_max)
    return v_max, v_min

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
    if R_min == 0:
        R_min = 1e-10
    v_R_min = (v_rot_univ(R_min) * r_sun / R_min - v_sun) * sin(l) * cos(b)
    v_R_max = (v_rot_univ(R_max) * r_sun / R_max - v_sun) * sin(l) * cos(b)
    v_max = np.maximum(v_R_min, v_R_max)
    v_min = np.minimum(v_R_min, v_R_max)
    return v_max, v_min

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
    v_R_min = (v_rot_linear(R_min) * r_sun / R_min - v_sun) * sin(l) * cos(b)
    v_R_max = (v_rot_linear(R_max) * r_sun / R_max - v_sun) * sin(l) * cos(b)
    v_max = np.maximum(v_R_min, v_R_max)
    v_min = np.minimum(v_R_min, v_R_max)
    return v_max, v_min

# Calculate V_max and V_min for given (l, b, h, r_gal)
# Polynomial rotation curve model
def calc_v_max_min_poly(l, b, h, r_gal):
    r_sun = 8.5
    v_sun = v_rot_poly(r_sun)
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
    v_R_min = (v_rot_poly(R_min) * r_sun / R_min - v_sun) * sin(l) * cos(b)
    v_R_max = (v_rot_poly(R_max) * r_sun / R_max - v_sun) * sin(l) * cos(b)
    v_max = np.maximum(v_R_min, v_R_max)
    v_min = np.minimum(v_R_min, v_R_max)
    return v_max, v_min

# Deviation Velocity
def calc_v_dev(
    l, b, h=5, r_gal=20, r_sun=8.5, v_sun=220, r_cut=0.5, model="univ", v_dev=0
):
    if model == "simple":
        print("model == simple")
        v_max, v_min = calc_v_max_min_simple(l, b, h, r_gal, r_sun, v_sun, r_cut)
    elif model == "univ":
        print("model == univ")
        v_max, v_min = calc_v_max_min_univ(l, b, h, r_gal)
    elif model == "linear":
        print("model == linear")
        v_max, v_min = calc_v_max_min_linear(l, b, h, r_gal)
    elif model == "poly":
        print("model == poly")
        v_max, v_min = calc_v_max_min_poly(l, b, h, r_gal)
    else:
        print("Invalid model")
        raise ValueError("Invalid model. Choose from 'simple', 'univ', 'linear', 'poly'.")
    return v_max + v_dev, v_min - v_dev