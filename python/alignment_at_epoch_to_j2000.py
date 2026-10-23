"""AlignmentAtEpoch -> J2000 的计算接口。"""
from pathlib import Path
import numpy as np
from time_utils import valid_epoch
from read_alignment_eop import read_alignment_eop
from fk5_j2000_to_fixed import fk5_j2000_to_fixed
from orbital_elements import elements_to_alignment


def alignment_at_epoch_to_j2000(coe, epoch, eop_file=None, *,
                                coordinate_epoch=None, mu=398600.4415,
                                use_polar_motion=True, use_nutation_corrections=False):
    """返回 (rJ, vJ, info)，位置 km、速度 km/s，数组形状均为 (3, 1)。

    coe=[a(km), e, i(deg), RAAN(deg), w(deg), M(deg)]；epoch 为 UTC 六元组。
    坐标系冻结在 coordinate_epoch，默认等于轨道历元；不进行轨道传播。
    保持 MATLAB 版本约定：速度不加入地球自转叉乘项。
    """
    epoch = valid_epoch(epoch)
    ce = valid_epoch(epoch if coordinate_epoch is None else coordinate_epoch)
    if not np.isscalar(mu) or not np.isreal(mu) or not np.isfinite(mu) or mu <= 0:
        raise ValueError('mu 必须是正的有限实数')
    for flag in (use_polar_motion, use_nutation_corrections):
        if not np.isscalar(flag) or flag not in (0, 1):
            raise ValueError('极移/章动修正开关必须为 True 或 False')
    r_a, v_a = elements_to_alignment(coe, mu)
    if eop_file is None:
        eop_file = Path(__file__).resolve().parent / 'data' / 'EOP-Last5Years.txt'
    q = read_alignment_eop(eop_file, ce)
    as2r = np.pi / 648000
    pm = np.array([q['xp'], q['yp']]) * as2r * use_polar_motion
    dn = np.array([q['dPsi'], q['dEpsilon']]) * as2r * use_nutation_corrections
    c, frame = fk5_j2000_to_fixed(ce, q['deltaAT'], q['deltaUT1'], pm, dn)
    options = dict(CoordinateEpoch=ce, Mu=mu, UsePolarMotion=bool(use_polar_motion),
                   UseNutationCorrections=bool(use_nutation_corrections))
    info = dict(OrbitEpochUTC=epoch, CoordinateEpochUTC=ce, Mu=mu, Model='IAU-76/FK5',
                EOP=q, rAlignment=r_a, vAlignment=v_a, C_J2000_to_Fixed=c,
                C_Fixed_to_J2000=c.T, Options=options, Frame=frame)
    return c.T @ r_a, c.T @ v_a, info
