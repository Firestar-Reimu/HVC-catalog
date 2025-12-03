/*
process_cube_c.c

C extension to speed up masking of a FITS data cube according to per-pixel
deviation-velocity ranges computed from Galactic coordinates (l,b).

This module exposes one function:

    process_mask_cube_numpy(data, vel_axis, l_map, b_map, model, n_sample, vdev)

Inputs (Python / numpy objects):
- data: 3D numpy array, shape (nz, ny, nx), dtype float64 (spectral axis first)
- vel_axis: 1D numpy array, length nz, dtype float64, velocities in km/s
- l_map: 2D numpy array, shape (ny, nx), dtype float64, Galactic longitude degrees
- b_map: 2D numpy array, shape (ny, nx), dtype float64, Galactic latitude degrees
- model: str, one of "simple","univ","linear","linear2","flat","power"
- n_sample: int (used by 'univ' and sampling in algorithms where needed)
- vdev: double, extra dev to add/subtract to v_max/v_min when computing ranges

Outputs (tuple):
(masked_cube, cube_gt_minvmax, cube_lt_maxvmin, vmax_map, vmin_map, global_min_vmax, global_max_vmin)

Notes:
- The function implements the rotation-curve and R_min/R_max logic internally (copied/ported
  from the earlier C/Python implementations) to avoid costly cross-language calls per-pixel.
- The function expects that vel_axis and spatial maps are prepared in Python (we do NOT call WCS here).
- All arrays are treated as double precision; other dtypes will be cast.
- The function modifies no input arrays; it returns new arrays.
*/

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#include <float.h>
#include <string.h>

/* Degree-based trig helpers */
static inline double deg_to_rad(double deg) { return deg * M_PI / 180.0; }
static inline double deg_sin(double deg) { return sin(deg_to_rad(deg)); }
static inline double deg_cos(double deg) { return cos(deg_to_rad(deg)); }
static inline double deg_tan(double deg) { return tan(deg_to_rad(deg)); }

/* Rotation models and R_min/R_max logic (translated from Python code) */

/* v_rot_simple (Wakker-like, parameterized by v_sun) */
static double v_rot_simple_model(double r, double v_sun_param) {
    if (r > 0.5) return v_sun_param;
    return 2.0 * v_sun_param * r;
}

/* v_rot_univ: universal curve with fixed internal r_sun = 8.15 (per user's requirement) */
static double v_rot_univ_model(double r) {
    const double a2 = 0.96;
    const double a3 = 1.62;
    const double r_sun_local = 8.15;

    double lambda_ = pow((a3 / 1.5), 5.0);
    double r_opt = a2 * r_sun_local;
    double x = r / r_opt;
    double log_lam = log10(lambda_);

    double term1 = 200.0 * pow(lambda_, 0.41);

    double top = 0.75 * exp(-0.4 * lambda_);
    double bot = 0.47 + 2.25 * pow(lambda_, 0.4);
    double term2 = sqrt(0.80 + 0.49 * log_lam + (top / bot));

    top = 1.97 * pow(x, 1.22);
    bot = pow(x*x + 0.61, 1.43);
    double term3 = (0.72 + 0.44 * log_lam) * (top / bot);

    top = x*x;
    bot = x*x + 2.25 * pow(lambda_, 0.4);
    double term4 = 1.6 * exp(-0.4 * lambda_) * (top / bot);

    double Vr = (term1 / term2) * sqrt(term3 + term4);
    return Vr;
}

/* v_rot_linear: uses r_sun_local = 8.122 */
static double v_rot_linear_model(double r) {
    const double r_sun_local = 8.122;
    return 229.0 - 1.7 * (r - r_sun_local);
}

/* v_rot_linear2: piecewise */
static double v_rot_linear2_model(double r) {
    if (r > 5.0) return 229.0 - 1.7 * r;
    return 44.1 * r;
}

/* v_rot_flat: constant v_sun */
static double v_rot_flat_model(double r, double v_sun_param) {
    (void)r;
    return v_sun_param;
}

/* v_rot_power: Russeil 2017 */
static double v_rot_power_model(double r) {
    const double r_sun_local = 8.34;
    const double v_sun_local = 240.0;
    return v_sun_local * 1.022 * pow(r / r_sun_local, 0.0803);
}

/* R_min/R_max calculators (1..4). These match the Python branching logic. */
static void calc_R_min_max_1(double l, double b, double h, double r_gal, double r_sun, double *Rmin, double *Rmax) {
    double L = fmod(l, 360.0);
    if (L < 0) L += 360.0;
    if (b == 0.0) b = 1e-10;
    double r_b = h / fabs(deg_tan(b));
    double r_l = sqrt(r_sun*r_sun + r_b*r_b - 2.0*r_sun*r_b*deg_cos(L));
    double cl = deg_cos(L);
    double y_0;
    if (r_b >= r_gal + r_sun) {
        if (cl >= 0) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_gal; }
        else { *Rmin = r_sun; *Rmax = r_gal; }
        return;
    } else if (2.0 * r_sun <= r_b && r_b < r_gal + r_sun) {
        y_0 = (r_gal*r_gal + r_sun*r_sun - r_b*r_b) / (2.0 * r_sun);
        if (cl >= (r_sun - y_0) / r_b) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_l; }
        else if (0.0 <= cl && cl < (r_sun - y_0) / r_b) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_gal; }
        else { *Rmin = r_sun; *Rmax = r_gal; }
        return;
    } else if (sqrt(r_gal*r_gal - r_sun*r_sun) <= r_b && r_b < 2.0 * r_sun) {
        y_0 = (r_gal*r_gal + r_sun*r_sun - r_b*r_b) / (2.0 * r_sun);
        if (cl >= r_b / (2.0 * r_sun)) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_sun; }
        else if ((r_sun - y_0) / r_b <= cl && cl < r_b / (2.0 * r_sun)) { *Rmin = r_sun; *Rmax = r_l; }
        else { *Rmin = r_sun; *Rmax = r_gal; }
        return;
    } else if (r_sun <= r_b && r_b < sqrt(r_gal*r_gal - r_sun*r_sun)) {
        y_0 = (r_gal*r_gal + r_sun*r_sun - r_b*r_b) / (2.0 * r_sun);
        if (cl >= r_b / (2.0 * r_sun)) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_sun; }
        else if (0.0 <= cl && cl < r_b / (2.0 * r_sun)) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_l; }
        else if ((r_sun - y_0) / r_b <= cl && cl < 0.0) { *Rmin = r_sun; *Rmax = r_l; }
        else { *Rmin = r_sun; *Rmax = r_gal; }
        return;
    } else if (r_gal - r_sun <= r_b && r_b < r_sun) {
        y_0 = (r_gal*r_gal + r_sun*r_sun - r_b*r_b) / (2.0 * r_sun);
        if (cl >= r_b / (2.0 * r_sun)) { *Rmin = r_l; *Rmax = r_sun; }
        else if (r_b / (2.0 * r_sun) <= cl && cl < r_b / r_sun) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_sun; }
        else if (0.0 <= cl && cl < r_b / (2.0 * r_sun)) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_l; }
        else if ((r_sun - y_0) / r_b <= cl && cl < 0.0) { *Rmin = r_sun; *Rmax = r_l; }
        else { *Rmin = r_sun; *Rmax = r_l; }
        return;
    } else {
        if (cl >= r_b / r_sun) { *Rmin = r_l; *Rmax = r_sun; }
        else if (r_b / (2.0 * r_sun) <= cl && cl < r_b / r_sun) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_sun; }
        else if (0.0 <= cl && cl < r_b / (2.0 * r_sun)) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_l; }
        else { *Rmin = r_sun; *Rmax = r_l; }
        return;
    }
}

static void calc_R_min_max_2(double l, double b, double h, double r_gal, double r_sun, double *Rmin, double *Rmax) {
    double L = fmod(l, 360.0);
    if (L < 0) L += 360.0;
    if (b == 0.0) b = 1e-10;
    double r_b = h / fabs(deg_tan(b));
    double r_l = sqrt(r_sun*r_sun + r_b*r_b - 2.0*r_sun*r_b*deg_cos(L));
    double cl = deg_cos(L);
    double y_0;
    if (r_b >= r_gal + r_sun) {
        if (cl >= 0) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_gal; }
        else { *Rmin = r_sun; *Rmax = r_gal; }
        return;
    } else if (2.0 * r_sun <= r_b && r_b < r_gal + r_sun) {
        y_0 = (r_gal*r_gal + r_sun*r_sun - r_b*r_b) / (2.0 * r_sun);
        if (cl >= (r_sun - y_0) / r_b) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_l; }
        else if (0.0 <= cl && cl < (r_sun - y_0) / r_b) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_gal; }
        else { *Rmin = r_sun; *Rmax = r_gal; }
        return;
    } else if (sqrt(r_gal*r_gal - r_sun*r_sun) <= r_b && r_b < 2.0 * r_sun) {
        y_0 = (r_gal*r_gal + r_sun*r_sun - r_b*r_b) / (2.0 * r_sun);
        if (cl >= r_b / (2.0 * r_sun)) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_sun; }
        else if ((r_sun - y_0) / r_b <= cl && cl < r_b / (2.0 * r_sun)) { *Rmin = r_sun; *Rmax = r_l; }
        else { *Rmin = r_sun; *Rmax = r_gal; }
        return;
    } else if (r_gal - r_sun <= r_b && r_b < sqrt(r_gal*r_gal - r_sun*r_sun)) {
        y_0 = (r_gal*r_gal + r_sun*r_sun - r_b*r_b) / (2.0 * r_sun);
        if (cl >= r_b / (2.0 * r_sun)) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_sun; }
        else if (0.0 <= cl && cl < r_b / (2.0 * r_sun)) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_l; }
        else if ((r_sun - y_0) / r_b <= cl && cl < 0.0) { *Rmin = r_sun; *Rmax = r_l; }
        else { *Rmin = r_sun; *Rmax = r_gal; }
        return;
    } else if (r_sun <= r_b && r_b < r_gal - r_sun) {
        if (cl >= r_b / (2.0 * r_sun)) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_sun; }
        else if (0.0 <= cl && cl < r_b / (2.0 * r_sun)) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_l; }
        else { *Rmin = r_sun; *Rmax = r_l; }
        return;
    } else {
        if (cl >= r_b / r_sun) { *Rmin = r_l; *Rmax = r_sun; }
        else if (r_b / (2.0 * r_sun) <= cl && cl < r_b / r_sun) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_sun; }
        else if (0.0 <= cl && cl < r_b / (2.0 * r_sun)) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_l; }
        else { *Rmin = r_sun; *Rmax = r_l; }
        return;
    }
}

static void calc_R_min_max_3(double l, double b, double h, double r_gal, double r_sun, double *Rmin, double *Rmax) {
    double L = fmod(l, 360.0);
    if (L < 0) L += 360.0;
    if (b == 0.0) b = 1e-10;
    double r_b = h / fabs(deg_tan(b));
    double r_l = sqrt(r_sun*r_sun + r_b*r_b - 2.0*r_sun*r_b*deg_cos(L));
    double cl = deg_cos(L);
    double y_0;
    if (r_b >= r_gal + r_sun) {
        if (cl >= 0) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_gal; }
        else { *Rmin = r_sun; *Rmax = r_gal; }
        return;
    } else if (sqrt(r_gal*r_gal - r_sun*r_sun) <= r_b && r_b < r_gal + r_sun) {
        y_0 = (r_gal*r_gal + r_sun*r_sun - r_b*r_b) / (2.0 * r_sun);
        if (cl >= (r_sun - y_0) / r_b) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_l; }
        else if (0.0 <= cl && cl < (r_sun - y_0) / r_b) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_gal; }
        else { *Rmin = r_sun; *Rmax = r_gal; }
        return;
    } else if (2.0 * r_sun <= r_b && r_b < sqrt(r_gal*r_gal - r_sun*r_sun)) {
        y_0 = (r_gal*r_gal + r_sun*r_sun - r_b*r_b) / (2.0 * r_sun);
        if (cl >= 0.0) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_l; }
        else if ((r_sun - y_0) / r_b <= cl && cl < 0.0) { *Rmin = r_sun; *Rmax = r_l; }
        else { *Rmin = r_sun; *Rmax = r_gal; }
        return;
    } else if (r_gal - r_sun <= r_b && r_b < 2.0 * r_sun) {
        y_0 = (r_gal*r_gal + r_sun*r_sun - r_b*r_b) / (2.0 * r_sun);
        if (cl >= r_b / (2.0 * r_sun)) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_sun; }
        else if (0.0 <= cl && cl < r_b / (2.0 * r_sun)) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_l; }
        else if ((r_sun - y_0) / r_b <= cl && cl < 0.0) { *Rmin = r_sun; *Rmax = r_l; }
        else { *Rmin = r_sun; *Rmax = r_gal; }
        return;
    } else if (r_sun <= r_b && r_b < r_gal - r_sun) {
        if (cl >= r_b / (2.0 * r_sun)) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_sun; }
        else if (0.0 <= cl && cl < r_b / (2.0 * r_sun)) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_l; }
        else { *Rmin = r_sun; *Rmax = r_l; }
        return;
    } else {
        if (cl >= r_b / r_sun) { *Rmin = r_l; *Rmax = r_sun; }
        else if (r_b / (2.0 * r_sun) <= cl && cl < r_b / r_sun) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_sun; }
        else if (0.0 <= cl && cl < r_b / (2.0 * r_sun)) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_l; }
        else { *Rmin = r_sun; *Rmax = r_l; }
        return;
    }
}

static void calc_R_min_max_4(double l, double b, double h, double r_gal, double r_sun, double *Rmin, double *Rmax) {
    double L = fmod(l, 360.0);
    if (L < 0) L += 360.0;
    if (b == 0.0) b = 1e-10;
    double r_b = h / fabs(deg_tan(b));
    double r_l = sqrt(r_sun*r_sun + r_b*r_b - 2.0*r_sun*r_b*deg_cos(L));
    double cl = deg_cos(L);
    double y_0;
    if (r_b >= r_gal + r_sun) {
        if (cl >= 0) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_gal; }
        else { *Rmin = r_sun; *Rmax = r_gal; }
        return;
    } else if (sqrt(r_gal*r_gal - r_sun*r_sun) <= r_b && r_b < r_gal + r_sun) {
        y_0 = (r_gal*r_gal + r_sun*r_sun - r_b*r_b) / (2.0 * r_sun);
        if (cl >= (r_sun - y_0) / r_b) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_l; }
        else if (0.0 <= cl && cl < (r_sun - y_0) / r_b) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_gal; }
        else { *Rmin = r_sun; *Rmax = r_gal; }
        return;
    } else if (r_gal - r_sun <= r_b && r_b < sqrt(r_gal*r_gal - r_sun*r_sun)) {
        y_0 = (r_gal*r_gal + r_sun*r_sun - r_b*r_b) / (2.0 * r_sun);
        if (cl >= 0.0) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_l; }
        else if ((r_sun - y_0) / r_b <= cl && cl < 0.0) { *Rmin = r_sun; *Rmax = r_l; }
        else { *Rmin = r_sun; *Rmax = r_gal; }
        return;
    } else if (2.0 * r_sun <= r_b && r_b < r_gal - r_sun) {
        if (cl >= 0.0) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_sun; }
        else { *Rmin = r_sun; *Rmax = r_l; }
        return;
    } else if (r_sun <= r_b && r_b < 2.0 * r_sun) {
        if (cl >= r_b / (2.0 * r_sun)) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_sun; }
        else if (0.0 <= cl && cl < r_b / (2.0 * r_sun)) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_l; }
        else { *Rmin = r_sun; *Rmax = r_l; }
        return;
    } else {
        if (cl >= r_b / r_sun) { *Rmin = r_l; *Rmax = r_sun; }
        else if (r_b / (2.0 * r_sun) <= cl && cl < r_b / r_sun) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_sun; }
        else if (0.0 <= cl && cl < r_b / (2.0 * r_sun)) { *Rmin = r_sun * fabs(deg_sin(L)); *Rmax = r_l; }
        else { *Rmin = r_sun; *Rmax = r_l; }
        return;
    }
}

/* Core calc_v_dev logic (internal C) matching Python signature defaults.
   We provide a function that returns v_max and v_min for a single (l,b).
   Default arguments set to match Python user's latest version:
   calc_v_dev(l, b, h=5, r_gal=20, r_sun=8.5, v_sun=220, model='univ', v_dev=50)
   Note: for 'univ' and 'linear' models internal r_sun values are fixed as per user's request.
*/
static void calc_v_dev_internal(double l, double b,
                                double h, double r_gal, double r_sun,
                                double v_sun_param, int n_sample,
                                const char *model, double v_dev,
                                double *out_vmax, double *out_vmin)
{
    double R_min = 0.0, R_max = 0.0;
    if (r_sun < r_gal && r_gal < 2.0 * r_sun) {
        calc_R_min_max_1(l,b,h,r_gal,r_sun,&R_min,&R_max);
    } else if (2.0 * r_sun < r_gal && r_gal < sqrt(5.0) * r_sun) {
        calc_R_min_max_2(l,b,h,r_gal,r_sun,&R_min,&R_max);
    } else if (sqrt(5.0) * r_sun < r_gal && r_gal < 3.0 * r_sun) {
        calc_R_min_max_3(l,b,h,r_gal,r_sun,&R_min,&R_max);
    } else {
        calc_R_min_max_4(l,b,h,r_gal,r_sun,&R_min,&R_max);
    }
    if (R_min == 0.0) R_min = 1e-10;
    if (R_max == 0.0) R_max = 1e-10;

    double s_l = deg_sin(l);
    double c_b = deg_cos(b);

    double vmax = -DBL_MAX;
    double vmin = DBL_MAX;

    if (strcmp(model, "simple") == 0) {
        double v1 = (v_rot_simple_model(R_min, v_sun_param) * r_sun / R_min - v_sun_param) * s_l * c_b;
        double v2 = (v_rot_simple_model(R_max, v_sun_param) * r_sun / R_max - v_sun_param) * s_l * c_b;
        vmax = (v1 > v2 ? v1 : v2) + v_dev;
        vmin = (v1 < v2 ? v1 : v2) - v_dev;
    } else if (strcmp(model, "univ") == 0) {
        /* use internal v_rot_univ_model that uses r_sun = 8.15 internally */
        for (int idx = 0; idx <= n_sample; ++idx) {
            double R = R_min + (R_max - R_min) * ((double)idx / (double)n_sample);
            if (R <= 0.0) R = 1e-10;
            double vr = v_rot_univ_model(R);
            double v = (vr * r_sun / R - v_sun_param) * s_l * c_b;
            if (v > vmax) vmax = v;
            if (v < vmin) vmin = v;
        }
        vmax += v_dev;
        vmin -= v_dev;
    } else if (strcmp(model, "linear") == 0) {
        for (int idx = 0; idx <= n_sample; ++idx) {
            double R = R_min + (R_max - R_min) * ((double)idx / (double)n_sample);
            if (R <= 0.0) R = 1e-10;
            double vr = v_rot_linear_model(R);
            double v = (vr * r_sun / R - v_sun_param) * s_l * c_b;
            if (v > vmax) vmax = v;
            if (v < vmin) vmin = v;
        }
        vmax += v_dev;
        vmin -= v_dev;
    } else if (strcmp(model, "linear2") == 0) {
        for (int idx = 0; idx <= n_sample; ++idx) {
            double R = R_min + (R_max - R_min) * ((double)idx / (double)n_sample);
            if (R <= 0.0) R = 1e-10;
            double vr = v_rot_linear2_model(R);
            double v = (vr * r_sun / R - v_sun_param) * s_l * c_b;
            if (v > vmax) vmax = v;
            if (v < vmin) vmin = v;
        }
        vmax += v_dev;
        vmin -= v_dev;
    } else if (strcmp(model, "power") == 0) {
        for (int idx = 0; idx <= n_sample; ++idx) {
            double R = R_min + (R_max - R_min) * ((double)idx / (double)n_sample);
            if (R <= 0.0) R = 1e-10;
            double vr = v_rot_power_model(R);
            double v = (vr * r_sun / R - v_sun_param) * s_l * c_b;
            if (v > vmax) vmax = v;
            if (v < vmin) vmin = v;
        }
        vmax += v_dev;
        vmin -= v_dev;
    } else if (strcmp(model, "flat") == 0) {
        for (int idx = 0; idx <= n_sample; ++idx) {
            double R = R_min + (R_max - R_min) * ((double)idx / (double)n_sample);
            if (R <= 0.0) R = 1e-10;
            double vr = v_rot_flat_model(R, v_sun_param);
            double v = (vr * r_sun / R - v_sun_param) * s_l * c_b;
            if (v > vmax) vmax = v;
            if (v < vmin) vmin = v;
        }
        vmax += v_dev;
        vmin -= v_dev;
    } else {
        /* default: univ */
        for (int idx = 0; idx <= n_sample; ++idx) {
            double R = R_min + (R_max - R_min) * ((double)idx / (double)n_sample);
            if (R <= 0.0) R = 1e-10;
            double vr = v_rot_univ_model(R);
            double v = (vr * r_sun / R - v_sun_param) * s_l * c_b;
            if (v > vmax) vmax = v;
            if (v < vmin) vmin = v;
        }
        vmax += v_dev;
        vmin -= v_dev;
    }

    *out_vmax = vmax;
    *out_vmin = vmin;
}

/* Main exposed function: process_mask_cube_numpy
   Signature in Python:
     process_mask_cube_numpy(data, vel_axis, l_map, b_map, model, n_sample=1000, vdev=0.0)

   All arrays expected to be numpy (ndarray) of dtype float64; data shape (nz,ny,nx)
*/
static PyObject *process_mask_cube_numpy(PyObject *self, PyObject *args, PyObject *kwds) {
    PyObject *data_obj = NULL;
    PyObject *vel_obj = NULL;
    PyObject *lmap_obj = NULL;
    PyObject *bmap_obj = NULL;
    const char *model = "univ";
    int n_sample = 1000;
    double vdev = 0.0;

    static char *kwlist[] = {"data", "vel_axis", "l_map", "b_map", "model", "n_sample", "vdev", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOO|sid", kwlist,
                                     &data_obj, &vel_obj, &lmap_obj, &bmap_obj,
                                     &model, &n_sample, &vdev)) {
        return NULL;
    }

    /* Convert inputs to numpy arrays (double, contiguous) */
    PyArrayObject *data_arr = (PyArrayObject *)PyArray_FROM_OTF(data_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    PyArrayObject *vel_arr  = (PyArrayObject *)PyArray_FROM_OTF(vel_obj,  NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    PyArrayObject *lmap_arr = (PyArrayObject *)PyArray_FROM_OTF(lmap_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    PyArrayObject *bmap_arr = (PyArrayObject *)PyArray_FROM_OTF(bmap_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);

    if (!data_arr || !vel_arr || !lmap_arr || !bmap_arr) {
        Py_XDECREF(data_arr); Py_XDECREF(vel_arr); Py_XDECREF(lmap_arr); Py_XDECREF(bmap_arr);
        PyErr_SetString(PyExc_TypeError, "Failed to convert inputs to required numpy dtypes (float64).");
        return NULL;
    }

    int ndim = PyArray_NDIM(data_arr);
    if (ndim != 3) {
        Py_DECREF(data_arr); Py_DECREF(vel_arr); Py_DECREF(lmap_arr); Py_DECREF(bmap_arr);
        PyErr_SetString(PyExc_ValueError, "data must be a 3D array with shape (nz,ny,nx).");
        return NULL;
    }

    npy_intp *shape = PyArray_SHAPE(data_arr);
    npy_intp nz = shape[0];
    npy_intp ny = shape[1];
    npy_intp nx = shape[2];

    if (PyArray_NDIM(vel_arr) != 1 || PyArray_DIM(vel_arr,0) != nz) {
        Py_DECREF(data_arr); Py_DECREF(vel_arr); Py_DECREF(lmap_arr); Py_DECREF(bmap_arr);
        PyErr_SetString(PyExc_ValueError, "vel_axis must be 1D of length nz.");
        return NULL;
    }
    if (PyArray_NDIM(lmap_arr) != 2 || PyArray_DIM(lmap_arr,0) != ny || PyArray_DIM(lmap_arr,1) != nx) {
        Py_DECREF(data_arr); Py_DECREF(vel_arr); Py_DECREF(lmap_arr); Py_DECREF(bmap_arr);
        PyErr_SetString(PyExc_ValueError, "l_map must be 2D with shape (ny,nx).");
        return NULL;
    }
    if (PyArray_NDIM(bmap_arr) != 2 || PyArray_DIM(bmap_arr,0) != ny || PyArray_DIM(bmap_arr,1) != nx) {
        Py_DECREF(data_arr); Py_DECREF(vel_arr); Py_DECREF(lmap_arr); Py_DECREF(bmap_arr);
        PyErr_SetString(PyExc_ValueError, "b_map must be 2D with shape (ny,nx).");
        return NULL;
    }

    /* Create output arrays (copies) */
    PyArrayObject *masked = (PyArrayObject *)PyArray_NewCopy(data_arr, NPY_CORDER);
    PyArrayObject *cube_gt = (PyArrayObject *)PyArray_NewCopy(masked, NPY_CORDER);
    PyArrayObject *cube_lt = (PyArrayObject *)PyArray_NewCopy(masked, NPY_CORDER);

    /* vmax_map and vmin_map */
    npy_intp dims2[2]; dims2[0]=ny; dims2[1]=nx;
    PyArrayObject *vmax_map = (PyArrayObject *)PyArray_SimpleNew(2, dims2, NPY_DOUBLE);
    PyArrayObject *vmin_map = (PyArrayObject *)PyArray_SimpleNew(2, dims2, NPY_DOUBLE);

    if (!masked || !cube_gt || !cube_lt || !vmax_map || !vmin_map) {
        Py_XDECREF(data_arr); Py_XDECREF(vel_arr); Py_XDECREF(lmap_arr); Py_XDECREF(bmap_arr);
        Py_XDECREF(masked); Py_XDECREF(cube_gt); Py_XDECREF(cube_lt);
        Py_XDECREF(vmax_map); Py_XDECREF(vmin_map);
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate output arrays.");
        return NULL;
    }

    double *data_ptr = (double *)PyArray_DATA(masked);
    double *vel_ptr = (double *)PyArray_DATA(vel_arr);
    double *lptr = (double *)PyArray_DATA(lmap_arr);
    double *bptr = (double *)PyArray_DATA(bmap_arr);
    double *vmax_ptr = (double *)PyArray_DATA(vmax_map);
    double *vmin_ptr = (double *)PyArray_DATA(vmin_map);

    /* strides: assume C-order contiguous (nz,ny,nx) */
    npy_intp stride0 = PyArray_STRIDES(masked)[0] / sizeof(double); /* stride between spec indices */
    npy_intp stride1 = PyArray_STRIDES(masked)[1] / sizeof(double];
    npy_intp stride2 = PyArray_STRIDES(masked)[2] / sizeof(double];

    /* We'll access masked as data_ptr[(k*ny + j)*nx + i] if contiguous C-order */
    /* But to be safe, compute element index via indices and shape. */
    /* Compute mask per pixel and set values to 0.0 */

    /* Initialize vmax_map/vmin_map to NaN */
    for (npy_intp jj=0; jj<ny; ++jj) {
        for (npy_intp ii=0; ii<nx; ++ii) {
            vmax_ptr[jj*nx + ii] = NAN;
            vmin_ptr[jj*nx + ii] = NAN;
        }
    }

    /* Per-pixel processing */
    for (npy_intp jj=0; jj<ny; ++jj) {
        for (npy_intp ii=0; ii<nx; ++ii) {
            double lval = lptr[jj*nx + ii];
            double bval = bptr[jj*nx + ii];
            if (!isfinite(lval) || !isfinite(bval)) {
                vmax_ptr[jj*nx + ii] = NAN;
                vmin_ptr[jj*nx + ii] = NAN;
                continue;
            }
            double v_max, v_min;
            /* call internal calc_v_dev; use defaults matching Python: h=5, r_gal=20, r_sun=8.5, v_sun=220 */
            calc_v_dev_internal(lval, bval, 5.0, 20.0, 8.5, 220.0, n_sample, model, vdev, &v_max, &v_min);
            vmax_ptr[jj*nx + ii] = v_max;
            vmin_ptr[jj*nx + ii] = v_min;
            if (!isfinite(v_max) || !isfinite(v_min)) continue;
            double vlow = v_min < v_max ? v_min : v_max;
            double vhigh = v_min < v_max ? v_max : v_min;

            /* iterate over spectral channels and mask if vel in [vlow, vhigh] */
            for (npy_intp kk=0; kk<nz; ++kk) {
                double vchan = vel_ptr[kk];
                if (vchan >= vlow && vchan <= vhigh) {
                    /* set masked[kk,jj,ii] = 0.0, cube_gt/cube_lt will be derived later */
                    /* compute flat index for (kk,jj,ii) in masked C-contiguous array */
                    npy_intp index = (kk * ny + jj) * nx + ii;
                    data_ptr[index] = 0.0;
                }
            }
        }
    }

    /* compute global thresholds ignoring NaNs */
    double global_min_vmax = INFINITY;
    double global_max_vmin = -INFINITY;
    for (npy_intp jj=0; jj<ny; ++jj) {
        for (npy_intp ii=0; ii<nx; ++ii) {
            double vm = vmax_ptr[jj*nx + ii];
            double vn = vmin_ptr[jj*nx + ii];
            if (isfinite(vm) && vm < global_min_vmax) global_min_vmax = vm;
            if (isfinite(vn) && vn > global_max_vmin) global_max_vmin = vn;
        }
    }
    if (!isfinite(global_min_vmax)) global_min_vmax = NAN;
    if (!isfinite(global_max_vmin)) global_max_vmin = NAN;

    /* produce cube_gt and cube_lt from masked (masked already applied to 'masked') */
    double *cube_gt_ptr = (double *)PyArray_DATA(cube_gt);
    double *cube_lt_ptr = (double *)PyArray_DATA(cube_lt);
    /* copy masked into both outputs then zero channels according to global thresholds */
    /* masked data is already in 'masked' array (data_ptr) */
    npy_intp total_elems = nz * ny * nx;
    for (npy_intp idx=0; idx<total_elems; ++idx) {
        cube_gt_ptr[idx] = data_ptr[idx];
        cube_lt_ptr[idx] = data_ptr[idx];
    }

    if (isfinite(global_min_vmax)) {
        for (npy_intp kk=0; kk<nz; ++kk) {
            double vchan = vel_ptr[kk];
            if (vchan <= global_min_vmax) {
                /* zero this channel in cube_gt */
                for (npy_intp jj=0; jj<ny; ++jj) {
                    for (npy_intp ii=0; ii<nx; ++ii) {
                        npy_intp index = (kk * ny + jj) * nx + ii;
                        cube_gt_ptr[index] = 0.0;
                    }
                }
            }
        }
    }
    if (isfinite(global_max_vmin)) {
        for (npy_intp kk=0; kk<nz; ++kk) {
            double vchan = vel_ptr[kk];
            if (vchan >= global_max_vmin) {
                /* zero this channel in cube_lt */
                for (npy_intp jj=0; jj<ny; ++jj) {
                    for (npy_intp ii=0; ii<nx; ++ii) {
                        npy_intp index = (kk * ny + jj) * nx + ii;
                        cube_lt_ptr[index] = 0.0;
                    }
                }
            }
        }
    }

    /* Build return tuple */
    PyObject *ret_masked = PyArray_Return(masked);
    PyObject *ret_gt = PyArray_Return(cube_gt);
    PyObject *ret_lt = PyArray_Return(cube_lt);
    PyObject *ret_vmax = PyArray_Return(vmax_map);
    PyObject *ret_vmin = PyArray_Return(vmin_map);
    PyObject *ret_global_min = PyFloat_FromDouble(global_min_vmax);
    PyObject *ret_global_max = PyFloat_FromDouble(global_max_vmin);

    /* cleanup references to inputs */
    Py_DECREF(data_arr); Py_DECREF(vel_arr); Py_DECREF(lmap_arr); Py_DECREF(bmap_arr);

    PyObject *result = Py_BuildValue("NNNNNN", ret_masked, ret_gt, ret_lt, ret_vmax, ret_vmin, Py_BuildValue("dd", global_min_vmax, global_max_vmin));
    /* The above embeds a tuple of two doubles as the 6th element (for convenience).
       Alternatively could return separately; consumers should unpack accordingly. */

    return result;
}

/* Method table */
static PyMethodDef ProcessCubeMethods[] = {
    {"process_mask_cube_numpy", (PyCFunction)process_mask_cube_numpy, METH_VARARGS | METH_KEYWORDS,
     "process_mask_cube_numpy(data, vel_axis, l_map, b_map, model='univ', n_sample=1000, vdev=0.0)\n"
     "Return (masked_cube, cube_gt, cube_lt, vmax_map, vmin_map, (global_min_vmax, global_max_vmin))."},
    {NULL, NULL, 0, NULL}
};

/* Module definition */
static struct PyModuleDef processcubemodule = {
    PyModuleDef_HEAD_INIT,
    "process_cube_c",
    "C-accelerated processing of FITS cubes by per-pixel dev-velocity masking",
    -1,
    ProcessCubeMethods
};

PyMODINIT_FUNC PyInit_process_cube_c(void) {
    PyObject *m;
    m = PyModule_Create(&processcubemodule);
    if (m == NULL) return NULL;
    /* Initialize numpy */
    import_array();
    return m;
}