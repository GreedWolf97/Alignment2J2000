"""UTC 输入校验与公历儒略日计算。"""
from datetime import date
import numpy as np


def real_vector(value, size, name):
    a = np.asarray(value)
    if (a.dtype.kind not in 'iuf' or a.size != size or a.ndim > 2
            or (a.ndim == 2 and 1 not in a.shape)):
        raise ValueError(f'{name} 必须是包含 {size} 个实数的向量')
    a = a.astype(float).reshape(-1)
    if not np.all(np.isfinite(a)):
        raise ValueError(f'{name} 必须是有限实数')
    return a


def valid_epoch(value):
    """返回 UTC [年, 月, 日, 时, 分, 秒]，不接受闰秒 60。"""
    t = real_vector(value, 6, 'UTC')
    if (np.any(t[:5] != np.fix(t[:5])) or not 0 <= t[3] < 24
            or not 0 <= t[4] < 60 or not 0 <= t[5] < 60):
        raise ValueError('UTC 格式应为 [Y M D h m s]，0 <= s < 60')
    date(*(int(x) for x in t[:3]))
    return t


def julian_date_midnight(utc):
    """公历当日 UTC 0 时的儒略日；保留日内秒单独计算以避免精度损失。"""
    return date(*(int(x) for x in utc[:3])).toordinal() + 1721424.5
