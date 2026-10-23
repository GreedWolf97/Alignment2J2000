"""右手系主动旋转矩阵（角度单位 rad）。"""
from math import cos, sin
import numpy as np


def rx(x):
    c, s = cos(x), sin(x)
    return np.array([[1., 0., 0.], [0., c, -s], [0., s, c]])


def ry(x):
    c, s = cos(x), sin(x)
    return np.array([[c, 0., s], [0., 1., 0.], [-s, 0., c]])


def rz(x):
    c, s = cos(x), sin(x)
    return np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]])
