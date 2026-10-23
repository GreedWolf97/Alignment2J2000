"""椭圆轨道根数转换到 Alignment 坐标系位置、速度。"""
from math import pi, sin, cos, sqrt
import numpy as np
from rotations import rx, rz
from time_utils import real_vector


def elements_to_alignment(coe, mu):
    """coe=[a(km), e, i(deg), RAAN(deg), w(deg), M(deg)]。"""
    c = real_vector(coe, 6, 'coe')
    a, e = c[:2]
    if a <= 0 or not 0 <= e < 1 or not 0 <= c[2] <= 180:
        raise ValueError('轨道根数必须满足 a>0, 0<=e<1, 0<=i<=180')
    m = (c[5] * pi / 180 + pi) % (2*pi) - pi
    lo, hi = -pi, pi
    for _ in range(60):
        eccentric = (lo + hi) / 2
        if eccentric - e*sin(eccentric) > m:
            hi = eccentric
        else:
            lo = eccentric
    eccentric = (lo + hi) / 2
    b = sqrt(1 - e**2)
    rp = a * np.array([cos(eccentric) - e, b*sin(eccentric), 0.])
    vp = sqrt(mu / a) / (1 - e*cos(eccentric)) * np.array([-sin(eccentric), b*cos(eccentric), 0.])
    q = rz(c[3]*pi/180) @ rx(c[2]*pi/180) @ rz(c[4]*pi/180)
    return (q @ rp).reshape(3, 1), (q @ vp).reshape(3, 1)
