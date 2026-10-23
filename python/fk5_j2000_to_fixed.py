"""IAU 1976 岁差、完整 106 项 IAU 1980 章动及 IERS 1996 春分点方程。"""
from functools import lru_cache
from pathlib import Path
from math import pi, sin, cos
import numpy as np
from rotations import rx, ry, rz
from time_utils import valid_epoch, julian_date_midnight


@lru_cache(maxsize=1)
def _nutation_table():
    # 原始数据来源：CelesTrak/fundamentals-of-astrodynamics/datalib。
    tab = np.loadtxt(Path(__file__).resolve().parent / 'data' / 'nut80.dat')
    if tab.ndim != 2 or tab.shape[0] != 106 or tab.shape[1] < 9 or not np.all(np.isfinite(tab)):
        raise ValueError('需要完整的 106 项 IAU 1980 章动表')
    tab.setflags(write=False)
    return tab


def fk5_j2000_to_fixed(utc, dat, dut1, pm, dn):
    """返回 J2000 -> ITRF 矩阵及中间量；pm、dn 的单位为 rad。"""
    utc = valid_epoch(utc)
    jd0 = julian_date_midnight(utc)
    s = utc[3] * 3600 + utc[4] * 60 + utc[5]
    t = ((jd0 - 2451545) + (s + dat + 32.184) / 86400) / 36525
    tu = ((jd0 - 2451545) + (s + dut1) / 86400) / 36525
    as2r = pi / 648000
    zeta = (2306.2181*t + 0.30188*t**2 + 0.017998*t**3) * as2r
    theta = (2004.3109*t - 0.42665*t**2 - 0.041833*t**3) * as2r
    z = (2306.2181*t + 1.09468*t**2 + 0.018203*t**3) * as2r
    p = rz(z) @ ry(-theta) @ rz(zeta)
    args = np.array([
        134.96298139*3600 + 1717915922.633*t + 31.310*t**2 + 0.064*t**3,
        357.52772333*3600 + 129596581.224*t - 0.577*t**2 - 0.012*t**3,
        93.27191028*3600 + 1739527263.137*t - 13.257*t**2 + 0.011*t**3,
        297.85036306*3600 + 1602961601.328*t - 6.891*t**2 + 0.019*t**3,
        125.04452222*3600 - 6962890.539*t + 7.455*t**2 + 0.008*t**3,
    ]) * as2r
    tab = _nutation_table()
    phase = tab[:, :5] @ args
    dpsi = np.sum((tab[:, 5] + tab[:, 6]*t) * np.sin(phase)) * 1e-4 * as2r + dn[0]
    deps = np.sum((tab[:, 7] + tab[:, 8]*t) * np.cos(phase)) * 1e-4 * as2r + dn[1]
    eps = (84381.448 - 46.8150*t - 0.00059*t**2 + 0.001813*t**3) * as2r
    n = rx(eps + deps) @ rz(dpsi) @ rx(-eps)
    gmst = ((67310.54841 + (876600*3600 + 8640184.812866)*tu
             + 0.093104*tu**2 - 6.2e-6*tu**3) * pi / 43200) % (2*pi)
    eq = dpsi * cos(eps)
    if jd0 + (s + dut1) / 86400 > 2450449.5:
        eq += (0.00264*sin(args[4]) + 0.000063*sin(2*args[4])) * as2r
    gast = gmst + eq
    w = ry(pm[0]) @ rx(pm[1])
    c = w @ rz(-gast) @ n @ p
    return c, dict(T_TT=t, GMSTRad=gmst, GASTRad=gast, dPsiRad=dpsi,
                   dEpsilonRad=deps, P=p, N=n, W=w)
