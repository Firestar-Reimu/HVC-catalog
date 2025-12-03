from rotation_curve_c import calc_v_dev
from Galaxy_Rotation_Model.py import calc_v_dev_py

v_dev_max, v_dev_min = calc_v_dev(230, 0, model='univ')
print(v_dev_max, v_dev_min)

v_dev_max1, v_dev_min1 = calc_v_dev(286, 0, h=4.0, r_gal=26.0, model='simple')
print(v_dev_max1, v_dev_min1)

v_dev_max2, v_dev_min2 = calc_v_dev(70.0, 0, h=5.0, r_gal=25.0, r_sun=8.5, n_sample=500, model='linear', dev=40.0)
print(v_dev_max2, v_dev_min2)

v_dev_max3, v_dev_min3 = calc_v_dev(120.0, 0, h=5.0, r_gal=25.0, r_sun=8.5, n_sample=500, model='linear2', dev=40.0)
print(v_dev_max3, v_dev_min3)
