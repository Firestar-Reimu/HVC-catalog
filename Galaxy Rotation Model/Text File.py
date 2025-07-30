import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord

h = 4
r_gal = 26
r_sun = 8.5
v_sun = 220

sin = lambda degrees: np.sin(np.deg2rad(degrees))
cos = lambda degrees: np.cos(np.deg2rad(degrees))
tan = lambda degrees: np.tan(np.deg2rad(degrees))

def v_rot_simple(r):
    v_rot = np.where(r > 0.5, 220, 440 * r)
    return v_rot

def calc_v_max_min_simple_vec(l, b, h, r_gal, r_sun=8.5, v_sun=220):
    # 保证所有输入为array
    l = np.asarray(l) % 360
    b = np.asarray(b)
    h = np.asarray(h)
    r_gal = np.asarray(r_gal)
    r_sun = float(r_sun)

    # 角度相关
    sin = lambda x: np.sin(np.deg2rad(x))
    cos = lambda x: np.cos(np.deg2rad(x))
    tan = lambda x: np.tan(np.deg2rad(x))
    cosl = cos(l)
    sinl = sin(l)
    cosb = cos(b)
    abs_sinl = np.abs(sinl)
    b_safe = np.where(b == 0, 1e-10, b)
    tanb = tan(b_safe)

    r_b = h / np.abs(tanb)
    r_l = np.sqrt(r_sun**2 + r_b**2 - 2 * r_sun * r_b * cosl)

    # 四种r_gal区域
    conds_main = [
        (r_sun < r_gal) & (r_gal < 2 * r_sun),  # r_sun < r_gal < 2 r_sun
        (2 * r_sun < r_gal) & (r_gal < np.sqrt(5) * r_sun),
        (np.sqrt(5) * r_sun < r_gal) & (r_gal < 3 * r_sun),
        (r_gal >= 3 * r_sun)
    ]

    # 对于每个主case，定义对应的R_min, R_max表达式（这里需要填入每个calc_R_min_max_X的np.select版）
    def R_min_max_case1():
        # r_sun < r_gal < 2 r_sun
        y_0 = (r_gal**2 + r_sun**2 - r_b**2) / (2 * r_sun)
        conds = [
            r_b >= r_gal + r_sun,
            (2 * r_sun <= r_b) & (r_b < r_gal + r_sun),
            (np.sqrt(r_gal**2 - r_sun**2) <= r_b) & (r_b < 2 * r_sun),
            (r_sun <= r_b) & (r_b < np.sqrt(r_gal**2 - r_sun**2)),
            (r_gal - r_sun <= r_b) & (r_b < r_sun),
            r_b < r_gal - r_sun
        ]
        # 这里每个分支再np.select一次，嵌套条件详见你的原分支
        # 例如，第1分支
        conds_1 = [cosl >= 0]
        R_min1 = np.select(conds_1, [r_sun * abs_sinl], default=r_sun)
        R_max1 = np.select(conds_1, [r_gal], default=r_gal)

        # 其它分支同理，详见你的代码
        # 这里只写第一分支，其它请仿照你的原始条件补充
        R_min = np.select(conds, [
            R_min1,      # 分支1
            r_sun,       # 分支2（举例，需补全）
            r_sun,       # 分支3（举例，需补全）
            r_sun,       # 分支4（举例，需补全）
            r_sun,       # 分支5（举例，需补全）
            r_sun        # 分支6（举例，需补全）
        ])
        R_max = np.select(conds, [
            R_max1,      # 分支1
            r_gal,       # 分支2（举例，需补全）
            r_gal,       # 分支3（举例，需补全）
            r_gal,       # 分支4（举例，需补全）
            r_gal,       # 分支5（举例，需补全）
            r_gal        # 分支6（举例，需补全）
        ])
        return R_min, R_max

    # 其它主case请仿照case1定义，如
    def R_min_max_case2():
        # TODO: 按你的calc_R_min_max_2逻辑依次写出6个分支的np.select
        return R_min, R_max

    def R_min_max_case3():
        # TODO: 按你的calc_R_min_max_3逻辑依次写出6个分支的np.select
        return R_min, R_max

    def R_min_max_case4():
        # TODO: 按你的calc_R_min_max_4逻辑依次写出6个分支的np.select
        return R_min, R_max

    # 用np.select选择哪一个主case
    R_min_all = np.select(conds_main, [
        R_min_max_case1()[0],
        R_min_max_case2()[0],
        R_min_max_case3()[0],
        R_min_max_case4()[0]
    ])
    R_max_all = np.select(conds_main, [
        R_min_max_case1()[1],
        R_min_max_case2()[1],
        R_min_max_case3()[1],
        R_min_max_case4()[1]
    ])

    v_max = (v_rot_simple(R_min_all) * r_sun / R_min_all - v_sun) * sinl * cosb
    v_min = (v_rot_simple(R_max_all) * r_sun / R_max_all - v_sun) * sinl * cosb
    return v_max, v_min

# 用法举例
# l, b, h, r_gal 可以是等长的np.array
# v_max, v_min = calc_v_max_min_simple_vec(l, b, h, r_gal)


import numpy as np

def case1_rmin_rmax_flexible(l, b, h, r_gal, r_sun):
    """
    支持 l, b 为标量、1d或2d数组，自动 meshgrid 展开，输出 shape 为 (len(b), len(l))
    其它参数也会自动广播。
    """
    sin = lambda x: np.sin(np.deg2rad(x))
    cos = lambda x: np.cos(np.deg2rad(x))
    tan = lambda x: np.tan(np.deg2rad(x))

    # 标量转为ndarray
    l = np.atleast_1d(l)
    b = np.atleast_1d(b)
    h = np.asarray(h)
    r_gal = np.asarray(r_gal)
    r_sun = float(r_sun)

    # meshgrid，如果都1d且长度>1
    if l.ndim == 1 and b.ndim == 1 and (l.size > 1 and b.size > 1):
        L, B = np.meshgrid(l, b, indexing='xy')
    else:
        # 至少有一个是标量或len为1，自动广播
        L, B = np.broadcast_arrays(l, b)

    B_safe = np.where(B == 0, 1e-10, B)
    r_b = h / np.abs(tan(B_safe))
    cosl = cos(L)
    sinl = sin(L)
    abs_sinl = np.abs(sinl)
    r_l = np.sqrt(r_sun**2 + r_b**2 - 2 * r_sun * r_b * cosl)
    y_0 = (r_gal**2 + r_sun**2 - r_b**2) / (2 * r_sun)

    # 分支1: r_b >= r_gal + r_sun
    cond1 = (r_b >= r_gal + r_sun)
    cond1a = (cosl >= 0)
    R_min_1 = np.where(cond1a, r_sun * abs_sinl, r_sun)
    R_max_1 = np.where(cond1a, r_gal, r_gal)

    # 分支2: 2 * r_sun <= r_b < r_gal + r_sun
    cond2 = (2 * r_sun <= r_b) & (r_b < r_gal + r_sun)
    cond2a = (cosl >= (r_sun - y_0) / r_b)
    cond2b = (0 <= cosl) & (cosl < (r_sun - y_0) / r_b)
    R_min_2 = np.select([cond2a, cond2b], [r_sun * abs_sinl, r_sun * abs_sinl], default=r_sun)
    R_max_2 = np.select([cond2a, cond2b], [r_l, r_gal], default=r_gal)

    # 分支3: sqrt(r_gal**2 - r_sun**2) <= r_b < 2 * r_sun
    cond3 = (np.sqrt(r_gal**2 - r_sun**2) <= r_b) & (r_b < 2 * r_sun)
    cond3a = (cosl >= r_b / (2 * r_sun))
    cond3b = ((r_sun - y_0) / r_b <= cosl) & (cosl < r_b / (2 * r_sun))
    R_min_3 = np.select([cond3a, cond3b], [r_sun * abs_sinl, r_sun], default=r_sun)
    R_max_3 = np.select([cond3a, cond3b], [r_sun, r_l], default=r_gal)

    # 分支4: r_sun <= r_b < sqrt(r_gal**2 - r_sun**2)
    cond4 = (r_sun <= r_b) & (r_b < np.sqrt(r_gal**2 - r_sun**2))
    cond4a = (cosl >= r_b / (2 * r_sun))
    cond4b = (0 <= cosl) & (cosl < r_b / (2 * r_sun))
    cond4c = ((r_sun - y_0) / r_b <= cosl) & (cosl < 0)
    R_min_4 = np.select([cond4a, cond4b, cond4c], [r_sun * abs_sinl, r_sun * abs_sinl, r_sun], default=r_sun)
    R_max_4 = np.select([cond4a, cond4b, cond4c], [r_sun, r_l, r_l], default=r_gal)

    # 分支5: r_gal - r_sun <= r_b < r_sun
    cond5 = (r_gal - r_sun <= r_b) & (r_b < r_sun)
    cond5a = (cosl >= r_b / (2 * r_sun))
    cond5b = (r_b / (2 * r_sun) <= cosl) & (cosl < r_b / r_sun)
    cond5c = (0 <= cosl) & (cosl < r_b / (2 * r_sun))
    cond5d = ((r_sun - y_0) / r_b <= cosl) & (cosl < 0)
    R_min_5 = np.select([cond5a, cond5b, cond5c, cond5d], [r_l, r_sun * abs_sinl, r_sun * abs_sinl, r_sun], default=r_sun)
    R_max_5 = np.select([cond5a, cond5b, cond5c, cond5d], [r_sun, r_sun, r_l, r_l], default=r_l)

    # 分支6: r_b < r_gal - r_sun
    cond6 = (r_b < r_gal - r_sun)
    cond6a = (cosl >= r_b / r_sun)
    cond6b = (r_b / (2 * r_sun) <= cosl) & (cosl < r_b / r_sun)
    cond6c = (0 <= cosl) & (cosl < r_b / (2 * r_sun))
    R_min_6 = np.select([cond6a, cond6b, cond6c], [r_l, r_sun * abs_sinl, r_sun * abs_sinl], default=r_sun)
    R_max_6 = np.select([cond6a, cond6b, cond6c], [r_sun, r_sun, r_l], default=r_l)

    # 总合并
    R_min = np.select([cond1, cond2, cond3, cond4, cond5, cond6],
                      [R_min_1, R_min_2, R_min_3, R_min_4, R_min_5, R_min_6], default=np.nan)
    R_max = np.select([cond1, cond2, cond3, cond4, cond5, cond6],
                      [R_max_1, R_max_2, R_max_3, R_max_4, R_max_5, R_max_6], default=np.nan)
    return R_min, R_max