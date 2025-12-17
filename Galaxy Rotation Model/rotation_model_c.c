/* rotation_curve.c
 *
 * C translation of the provided Python "Galaxy_Rotation_Model.py".
 * I preserved the Python comments that immediately precede each calc_R_min_max_i()
 * and each v_rot_xxx() block (excluding the "# %%" separators), as you requested.
 *
 * Exposes:
 *   calc_v_dev(l, b, h=5, r_gal=20, r_sun=8.5, v_sun=220, model='univ', v_dev=50)
 *
 * Build with:
 *   python setup.py build_ext --inplace
 *
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>
#include <float.h>
#include <string.h>

/* Helpers: degree-based trig */
static inline double deg_to_rad(double deg) { return deg * M_PI / 180.0; }
static inline double deg_sin(double deg) { return sin(deg_to_rad(deg)); }
static inline double deg_cos(double deg) { return cos(deg_to_rad(deg)); }
static inline double deg_tan(double deg) { return tan(deg_to_rad(deg)); }

/* ------------------------------------------------------------------
# r_sun < r_gal < 2 r_sun
*/
static void calc_R_min_max_1(double l, double b, double h, double r_gal, double r_sun, double *Rmin, double *Rmax) {
    double L = fmod(l, 360.0);
    if (L < 0.0) L += 360.0;
    if (b == 0.0) b = 1e-10;
    double r_b = h / fabs(deg_tan(b));
    double r_l = sqrt(r_sun*r_sun + r_b*r_b - 2.0 * r_sun * r_b * deg_cos(L));
    double cl = deg_cos(L);
    double y_0;
    if (r_b >= r_gal + r_sun) { /* 1 */
        if (cl >= 0.0) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_gal;
        } else {
            *Rmin = r_sun;
            *Rmax = r_gal;
        }
    } else if (2.0 * r_sun <= r_b && r_b < r_gal + r_sun) { /* 2 */
        y_0 = (r_gal*r_gal + r_sun*r_sun - r_b*r_b) / (2.0 * r_sun);
        if (cl >= (r_sun - y_0) / r_b) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_l;
        } else if (0.0 <= cl && cl < (r_sun - y_0) / r_b) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_gal;
        } else {
            *Rmin = r_sun;
            *Rmax = r_gal;
        }
    } else if (sqrt(r_gal*r_gal - r_sun*r_sun) <= r_b && r_b < 2.0 * r_sun) { /* 3 */
        y_0 = (r_gal*r_gal + r_sun*r_sun - r_b*r_b) / (2.0 * r_sun);
        if (cl >= r_b / (2.0 * r_sun)) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_sun;
        } else if ((r_sun - y_0) / r_b <= cl && cl < r_b / (2.0 * r_sun)) {
            *Rmin = r_sun;
            *Rmax = r_l;
        } else {
            *Rmin = r_sun;
            *Rmax = r_gal;
        }
    } else if (r_sun <= r_b && r_b < sqrt(r_gal*r_gal - r_sun*r_sun)) { /* 4 */
        y_0 = (r_gal*r_gal + r_sun*r_sun - r_b*r_b) / (2.0 * r_sun);
        if (cl >= r_b / (2.0 * r_sun)) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_sun;
        } else if (0.0 <= cl && cl < r_b / (2.0 * r_sun)) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_l;
        } else if ((r_sun - y_0) / r_b <= cl && cl < 0.0) {
            *Rmin = r_sun;
            *Rmax = r_l;
        } else {
            *Rmin = r_sun;
            *Rmax = r_gal;
        }
    } else if (r_gal - r_sun <= r_b && r_b < r_sun) { /* 5 */
        y_0 = (r_gal*r_gal + r_sun*r_sun - r_b*r_b) / (2.0 * r_sun);
        if (cl >= r_b / (2.0 * r_sun)) {
            *Rmin = r_l;
            *Rmax = r_sun;
        } else if (r_b / (2.0 * r_sun) <= cl && cl < r_b / r_sun) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_sun;
        } else if (0.0 <= cl && cl < r_b / (2.0 * r_sun)) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_l;
        } else if ((r_sun - y_0) / r_b <= cl && cl < 0.0) {
            *Rmin = r_sun;
            *Rmax = r_l;
        } else {
            *Rmin = r_sun;
            *Rmax = r_l;
        }
    } else { /* 6 */
        if (cl >= r_b / r_sun) {
            *Rmin = r_l;
            *Rmax = r_sun;
        } else if (r_b / (2.0 * r_sun) <= cl && cl < r_b / r_sun) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_sun;
        } else if (0.0 <= cl && cl < r_b / (2.0 * r_sun)) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_l;
        } else {
            *Rmin = r_sun;
            *Rmax = r_l;
        }
    }
}

/* ------------------------------------------------------------------
# 2 r_sun < r_gal < sqrt(5) r_sun
*/
static void calc_R_min_max_2(double l, double b, double h, double r_gal, double r_sun, double *Rmin, double *Rmax) {
    double L = fmod(l, 360.0);
    if (L < 0.0) L += 360.0;
    if (b == 0.0) b = 1e-10;
    double r_b = h / fabs(deg_tan(b));
    double r_l = sqrt(r_sun*r_sun + r_b*r_b - 2.0 * r_sun * r_b * deg_cos(L));
    double cl = deg_cos(L);
    double y_0;
    if (r_b >= r_gal + r_sun) { /* 1 */
        if (cl >= 0.0) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_gal;
        } else {
            *Rmin = r_sun;
            *Rmax = r_gal;
        }
    } else if (2.0 * r_sun <= r_b && r_b < r_gal + r_sun) { /* 2 */
        y_0 = (r_gal*r_gal + r_sun*r_sun - r_b*r_b) / (2.0 * r_sun);
        if (cl >= (r_sun - y_0) / r_b) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_l;
        } else if (0.0 <= cl && cl < (r_sun - y_0) / r_b) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_gal;
        } else {
            *Rmin = r_sun;
            *Rmax = r_gal;
        }
    } else if (sqrt(r_gal*r_gal - r_sun*r_sun) <= r_b && r_b < 2.0 * r_sun) { /* 3 */
        y_0 = (r_gal*r_gal + r_sun*r_sun - r_b*r_b) / (2.0 * r_sun);
        if (cl >= r_b / (2.0 * r_sun)) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_sun;
        } else if ((r_sun - y_0) / r_b <= cl && cl < r_b / (2.0 * r_sun)) {
            *Rmin = r_sun;
            *Rmax = r_l;
        } else {
            *Rmin = r_sun;
            *Rmax = r_gal;
        }
    } else if (r_gal - r_sun <= r_b && r_b < sqrt(r_gal*r_gal - r_sun*r_sun)) { /* 4 */
        y_0 = (r_gal*r_gal + r_sun*r_sun - r_b*r_b) / (2.0 * r_sun);
        if (cl >= r_b / (2.0 * r_sun)) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_sun;
        } else if (0.0 <= cl && cl < r_b / (2.0 * r_sun)) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_l;
        } else if ((r_sun - y_0) / r_b <= cl && cl < 0.0) {
            *Rmin = r_sun;
            *Rmax = r_l;
        } else {
            *Rmin = r_sun;
            *Rmax = r_gal;
        }
    } else if (r_sun <= r_b && r_b < r_gal - r_sun) { /* 5 */
        if (cl >= r_b / (2.0 * r_sun)) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_sun;
        } else if (0.0 <= cl && cl < r_b / (2.0 * r_sun)) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_l;
        } else {
            *Rmin = r_sun;
            *Rmax = r_l;
        }
    } else { /* 6 */
        if (cl >= r_b / r_sun) {
            *Rmin = r_l;
            *Rmax = r_sun;
        } else if (r_b / (2.0 * r_sun) <= cl && cl < r_b / r_sun) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_sun;
        } else if (0.0 <= cl && cl < r_b / (2.0 * r_sun)) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_l;
        } else {
            *Rmin = r_sun;
            *Rmax = r_l;
        }
    }
}

/* ------------------------------------------------------------------
# sqrt(5) r_sun < r_gal < 3 r_sun
*/
static void calc_R_min_max_3(double l, double b, double h, double r_gal, double r_sun, double *Rmin, double *Rmax) {
    double L = fmod(l, 360.0);
    if (L < 0.0) L += 360.0;
    if (b == 0.0) b = 1e-10;
    double r_b = h / fabs(deg_tan(b));
    double r_l = sqrt(r_sun*r_sun + r_b*r_b - 2.0 * r_sun * r_b * deg_cos(L));
    double cl = deg_cos(L);
    double y_0;
    if (r_b >= r_gal + r_sun) { /* 1 */
        if (cl >= 0.0) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_gal;
        } else {
            *Rmin = r_sun;
            *Rmax = r_gal;
        }
    } else if (sqrt(r_gal*r_gal - r_sun*r_sun) <= r_b && r_b < r_gal + r_sun) { /* 2 */
        y_0 = (r_gal*r_gal + r_sun*r_sun - r_b*r_b) / (2.0 * r_sun);
        if (cl >= (r_sun - y_0) / r_b) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_l;
        } else if (0.0 <= cl && cl < (r_sun - y_0) / r_b) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_gal;
        } else {
            *Rmin = r_sun;
            *Rmax = r_gal;
        }
    } else if (2.0 * r_sun <= r_b && r_b < sqrt(r_gal*r_gal - r_sun*r_sun)) { /* 3 */
        y_0 = (r_gal*r_gal + r_sun*r_sun - r_b*r_b) / (2.0 * r_sun);
        if (cl >= 0.0) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_l;
        } else if ((r_sun - y_0) / r_b <= cl && cl < 0.0) {
            *Rmin = r_sun;
            *Rmax = r_l;
        } else {
            *Rmin = r_sun;
            *Rmax = r_gal;
        }
    } else if (r_gal - r_sun <= r_b && r_b < 2.0 * r_sun) { /* 4 */
        y_0 = (r_gal*r_gal + r_sun*r_sun - r_b*r_b) / (2.0 * r_sun);
        if (cl >= r_b / (2.0 * r_sun)) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_sun;
        } else if (0.0 <= cl && cl < r_b / (2.0 * r_sun)) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_l;
        } else if ((r_sun - y_0) / r_b <= cl && cl < 0.0) {
            *Rmin = r_sun;
            *Rmax = r_l;
        } else {
            *Rmin = r_sun;
            *Rmax = r_gal;
        }
    } else if (r_sun <= r_b && r_b < r_gal - r_sun) { /* 5 */
        if (cl >= r_b / (2.0 * r_sun)) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_sun;
        } else if (0.0 <= cl && cl < r_b / (2.0 * r_sun)) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_l;
        } else {
            *Rmin = r_sun;
            *Rmax = r_l;
        }
    } else { /* 6 */
        if (cl >= r_b / r_sun) {
            *Rmin = r_l;
            *Rmax = r_sun;
        } else if (r_b / (2.0 * r_sun) <= cl && cl < r_b / r_sun) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_sun;
        } else if (0.0 <= cl && cl < r_b / (2.0 * r_sun)) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_l;
        } else {
            *Rmin = r_sun;
            *Rmax = r_l;
        }
    }
}

/* ------------------------------------------------------------------
# r_gal >= 3 r_sun
*/
static void calc_R_min_max_4(double l, double b, double h, double r_gal, double r_sun, double *Rmin, double *Rmax) {
    double L = fmod(l, 360.0);
    if (L < 0.0) L += 360.0;
    if (b == 0.0) b = 1e-10;
    double r_b = h / fabs(deg_tan(b));
    double r_l = sqrt(r_sun*r_sun + r_b*r_b - 2.0 * r_sun * r_b * deg_cos(L));
    double cl = deg_cos(L);
    double y_0;
    if (r_b >= r_gal + r_sun) { /* 1 */
        if (cl >= 0.0) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_gal;
        } else {
            *Rmin = r_sun;
            *Rmax = r_gal;
        }
    } else if (sqrt(r_gal*r_gal - r_sun*r_sun) <= r_b && r_b < r_gal + r_sun) { /* 2 */
        y_0 = (r_gal*r_gal + r_sun*r_sun - r_b*r_b) / (2.0 * r_sun);
        if (cl >= (r_sun - y_0) / r_b) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_l;
        } else if (0.0 <= cl && cl < (r_sun - y_0) / r_b) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_gal;
        } else {
            *Rmin = r_sun;
            *Rmax = r_gal;
        }
    } else if (r_gal - r_sun <= r_b && r_b < sqrt(r_gal*r_gal - r_sun*r_sun)) { /* 3 */
        y_0 = (r_gal*r_gal + r_sun*r_sun - r_b*r_b) / (2.0 * r_sun);
        if (cl >= 0.0) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_l;
        } else if ((r_sun - y_0) / r_b <= cl && cl < 0.0) {
            *Rmin = r_sun;
            *Rmax = r_l;
        } else {
            *Rmin = r_sun;
            *Rmax = r_gal;
        }
    } else if (2.0 * r_sun <= r_b && r_b < r_gal - r_sun) { /* 4 */
        if (cl >= 0.0) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_sun;
        } else {
            *Rmin = r_sun;
            *Rmax = r_l;
        }
    } else if (r_sun <= r_b && r_b < 2.0 * r_sun) { /* 5 */
        if (cl >= r_b / (2.0 * r_sun)) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_sun;
        } else if (0.0 <= cl && cl < r_b / (2.0 * r_sun)) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_l;
        } else {
            *Rmin = r_sun;
            *Rmax = r_l;
        }
    } else { /* 6 */
        if (cl >= r_b / r_sun) {
            *Rmin = r_l;
            *Rmax = r_sun;
        } else if (r_b / (2.0 * r_sun) <= cl && cl < r_b / r_sun) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_sun;
        } else if (0.0 <= cl && cl < r_b / (2.0 * r_sun)) {
            *Rmin = r_sun * fabs(deg_sin(L));
            *Rmax = r_l;
        } else {
            *Rmin = r_sun;
            *Rmax = r_l;
        }
    }
}

/* ------------------------------------------------------------------
# Simple rotation curve model
# V_rot = 220 km/s where R > 0.5 kpc, else V_rot = 440 * R
# Wakker 1991: 1991A&A...250..499W
*/
static double v_rot_simple(double r, double v_sun_param, double r_cut_param) {
    if (r > r_cut_param) return v_sun_param;
    return (1.0 / r_cut_param) * v_sun_param * r;
}

/* ------------------------------------------------------------------
# Universal rotation curve model
# Persic et al. 1996: 1996MNRAS.281...27P
# Reid et al. 2019: 2019ApJ...885..131R
*/
static double v_rot_univ(double r) {
    /* Python sets r_sun = 8.15 inside this function (fixed). Follow that. */
    const double r_sun_local = 8.15;
    const double a2 = 0.96;
    const double a3 = 1.62;

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

/* ------------------------------------------------------------------
# Linear rotation curve model
# V_rot = 229 - 1.7 * (R - R_sun)
# Eilers et al. 2019: 2019ApJ...871...120E
*/
static double v_rot_linear(double r) {
    const double r_sun_local = 8.122;
    return 229.0 - 1.7 * (r - r_sun_local);
}

/* ------------------------------------------------------------------
# Power law rotation curve model
# V_rot = 240 * 1.022 * (R / 8.34)^0.0803
# Russeil et al. 2017: 2017A&A...601L...5R
*/
static double v_rot_power(double r) {
    const double r_sun_local = 8.34;
    const double v_sun_local = 240.0;
    return v_sun_local * 1.022 * pow(r / r_sun_local, 0.0803);
}

/* ------------------------------------------------------------------
   calc_v_max_min_simple: mirrors Python calc_v_max_min_simple
*/
static void calc_v_max_min_simple(double l, double b, double h, double r_gal, double r_sun_param, double v_sun_param, double r_cut_param, double *vmax, double *vmin) {
    double v_sun_local = v_rot_simple(r_sun_param, v_sun_param, r_cut_param);
    double Rmin = 0.0, Rmax = 0.0;
    if (r_sun_param < r_gal && r_gal < 2.0 * r_sun_param) {
        calc_R_min_max_1(l,b,h,r_gal,r_sun_param,&Rmin,&Rmax);
    } else if (2.0 * r_sun_param < r_gal && r_gal < sqrt(5.0) * r_sun_param) {
        calc_R_min_max_2(l,b,h,r_gal,r_sun_param,&Rmin,&Rmax);
    } else if (sqrt(5.0) * r_sun_param < r_gal && r_gal < 3.0 * r_sun_param) {
        calc_R_min_max_3(l,b,h,r_gal,r_sun_param,&Rmin,&Rmax);
    } else {
        calc_R_min_max_4(l,b,h,r_gal,r_sun_param,&Rmin,&Rmax);
    }
    double s_l = deg_sin(l);
    double c_b = deg_cos(b);
    double vRmin = (v_rot_simple(Rmin, v_sun_param, r_cut_param) * r_sun_param / Rmin - v_sun_local) * s_l * c_b;
    double vRmax = (v_rot_simple(Rmax, v_sun_param, r_cut_param) * r_sun_param / Rmax - v_sun_local) * s_l * c_b;
    *vmax = (vRmin > vRmax) ? vRmin : vRmax;
    *vmin = (vRmin < vRmax) ? vRmin : vRmax;
}

/* ------------------------------------------------------------------
   calc_v_max_min_univ: mirrors Python calc_v_max_min_univ
*/
static void calc_v_max_min_univ(double l, double b, double h, double r_gal, double *vmax, double *vmin) {
    const double r_sun_local = 8.15;
    double v_sun_local = v_rot_univ(r_sun_local);
    double Rmin = 0.0, Rmax = 0.0;
    if (r_sun_local < r_gal && r_gal < 2.0 * r_sun_local) {
        calc_R_min_max_1(l,b,h,r_gal,r_sun_local,&Rmin,&Rmax);
    } else if (2.0 * r_sun_local < r_gal && r_gal < sqrt(5.0) * r_sun_local) {
        calc_R_min_max_2(l,b,h,r_gal,r_sun_local,&Rmin,&Rmax);
    } else if (sqrt(5.0) * r_sun_local < r_gal && r_gal < 3.0 * r_sun_local) {
        calc_R_min_max_3(l,b,h,r_gal,r_sun_local,&Rmin,&Rmax);
    } else {
        calc_R_min_max_4(l,b,h,r_gal,r_sun_local,&Rmin,&Rmax);
    }
    double s_l = deg_sin(l);
    double c_b = deg_cos(b);
    double vRmin = (v_rot_univ(Rmin) * r_sun_local / Rmin - v_sun_local) * s_l * c_b;
    double vRmax = (v_rot_univ(Rmax) * r_sun_local / Rmax - v_sun_local) * s_l * c_b;
    *vmax = (vRmin > vRmax) ? vRmin : vRmax;
    *vmin = (vRmin < vRmax) ? vRmin : vRmax;
}

/* ------------------------------------------------------------------
   calc_v_max_min_linear: mirrors Python calc_v_max_min_linear
*/
static void calc_v_max_min_linear(double l, double b, double h, double r_gal, double *vmax, double *vmin) {
    const double r_sun_local = 8.122;
    double v_sun_local = v_rot_linear(r_sun_local);
    double Rmin = 0.0, Rmax = 0.0;
    if (r_sun_local < r_gal && r_gal < 2.0 * r_sun_local) {
        calc_R_min_max_1(l,b,h,r_gal,r_sun_local,&Rmin,&Rmax);
    } else if (2.0 * r_sun_local < r_gal && r_gal < sqrt(5.0) * r_sun_local) {
        calc_R_min_max_2(l,b,h,r_gal,r_sun_local,&Rmin,&Rmax);
    } else if (sqrt(5.0) * r_sun_local < r_gal && r_gal < 3.0 * r_sun_local) {
        calc_R_min_max_3(l,b,h,r_gal,r_sun_local,&Rmin,&Rmax);
    } else {
        calc_R_min_max_4(l,b,h,r_gal,r_sun_local,&Rmin,&Rmax);
    }
    if (Rmin == 0.0) Rmin = 1e-10;
    if (Rmax == 0.0) Rmax = 1e-10;
    double s_l = deg_sin(l);
    double c_b = deg_cos(b);
    double vRmin = (v_rot_linear(Rmin) * r_sun_local / Rmin - v_sun_local) * s_l * c_b;
    double vRmax = (v_rot_linear(Rmax) * r_sun_local / Rmax - v_sun_local) * s_l * c_b;
    *vmax = (vRmin > vRmax) ? vRmin : vRmax;
    *vmin = (vRmin < vRmax) ? vRmin : vRmax;
}

/* ------------------------------------------------------------------
   calc_v_max_min_power: mirrors Python calc_v_max_min_power
*/
static void calc_v_max_min_power(double l, double b, double h, double r_gal, double *vmax, double *vmin) {
    const double r_sun_local = 8.34;
    double v_sun_local = v_rot_power(r_sun_local);
    double Rmin = 0.0, Rmax = 0.0;
    if (r_sun_local < r_gal && r_gal < 2.0 * r_sun_local) {
        calc_R_min_max_1(l,b,h,r_gal,r_sun_local,&Rmin,&Rmax);
    } else if (2.0 * r_sun_local < r_gal && r_gal < sqrt(5.0) * r_sun_local) {
        calc_R_min_max_2(l,b,h,r_gal,r_sun_local,&Rmin,&Rmax);
    } else if (sqrt(5.0) * r_sun_local < r_gal && r_gal < 3.0 * r_sun_local) {
        calc_R_min_max_3(l,b,h,r_gal,r_sun_local,&Rmin,&Rmax);
    } else {
        calc_R_min_max_4(l,b,h,r_gal,r_sun_local,&Rmin,&Rmax);
    }
    double s_l = deg_sin(l);
    double c_b = deg_cos(b);
    double vRmin = (v_rot_power(Rmin) * r_sun_local / Rmin - v_sun_local) * s_l * c_b;
    double vRmax = (v_rot_power(Rmax) * r_sun_local / Rmax - v_sun_local) * s_l * c_b;
    *vmax = (vRmin > vRmax) ? vRmin : vRmax;
    *vmin = (vRmin < vRmax) ? vRmin : vRmax;
}

/* ------------------------------------------------------------------
   Public Python wrapper:
   calc_v_dev(l, b, h=5, r_gal=20, r_sun=8.5, v_sun=220, model='univ', v_dev=50)
*/
static PyObject* py_calc_v_dev(PyObject *self, PyObject *args, PyObject *kwds) {
    static char *kwlist[] = {"l", "b", "h", "r_gal", "r_sun", "v_sun", "model", "v_dev", NULL};

    double l, b;
    double h = 5.0;
    double r_gal = 20.0;
    double r_sun = 8.5;
    double v_sun_param = 220.0;
    double r_cut_param = 0.5;
    PyObject *model_obj = NULL;
    double v_dev = 50.0;

    /* parse: required l,b; optional h, r_gal, r_sun, v_sun, model, v_dev */
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "dd|ddddOd", kwlist,
                                     &l, &b, &h, &r_gal, &r_sun, &v_sun_param, &model_obj, &v_dev)) {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments for calc_v_dev");
        return NULL;
    }

    const char *model = "univ";
    if (model_obj && model_obj != Py_None) {
        model = PyUnicode_AsUTF8(PyObject_Str(model_obj));
        if (!model) {
            PyErr_SetString(PyExc_TypeError, "model must be a string");
            return NULL;
        }
    }

    double v_max = 0.0, v_min = 0.0;
    if (strcmp(model, "simple") == 0) {
        calc_v_max_min_simple(l, b, h, r_gal, r_sun, v_sun_param, r_cut_param, &v_max, &v_min);
    } else if (strcmp(model, "univ") == 0) {
        calc_v_max_min_univ(l, b, h, r_gal, &v_max, &v_min);
    } else if (strcmp(model, "linear") == 0) {
        calc_v_max_min_linear(l, b, h, r_gal, &v_max, &v_min);
    } else if (strcmp(model, "power") == 0) {
        calc_v_max_min_power(l, b, h, r_gal, &v_max, &v_min);
    } else {
        PyErr_SetString(PyExc_ValueError, "Invalid model. Choose from 'simple', 'univ', 'linear', 'power'.");
        return NULL;
    }

    double out_max = v_max + v_dev;
    double out_min = v_min - v_dev;
    PyObject *ret = Py_BuildValue("dd", out_max, out_min);
    return ret;
}

/* Module method table */
static PyMethodDef RotationCurveMethods[] = {
    {"calc_v_dev", (PyCFunction)py_calc_v_dev, METH_VARARGS | METH_KEYWORDS,
     "calc_v_dev(l, b, h=5, r_gal=20, r_sun=8.5, v_sun=220, model='univ', v_dev=50)"},
    {NULL, NULL, 0, NULL}
};

/* Module definition */
static struct PyModuleDef rotationmodelmodule = {
    PyModuleDef_HEAD_INIT,
    "rotation_model_c",
    "C implementation of calc_v_dev translated from Galaxy_Rotation_Model.py",
    -1,
    RotationCurveMethods
};

PyMODINIT_FUNC PyInit_rotation_model_c(void) {
    return PyModule_Create(&rotationmodelmodule);
}