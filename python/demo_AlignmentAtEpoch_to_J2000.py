"""修改输入区后直接运行：python demo_AlignmentAtEpoch_to_J2000.py。"""
from pathlib import Path
import numpy as np
from alignment_at_epoch_to_j2000 import alignment_at_epoch_to_j2000


def main():
    # 输入区：a(m), e, i(deg), RAAN(deg), w(deg), M(deg)。
    coe_m = np.array([6922975, 0.001, 39.123, 354.453, 95.59, 63.36])
    # UTC [年, 月, 日, 时, 分, 秒]，秒可以含小数。
    epoch_utc = [2026, 10, 1, 4, 55, 22.190]
    coordinate_epoch_utc = epoch_utc
    eop_file = Path(__file__).resolve().parent / 'data' / 'EOP-Last5Years.txt'

    coe_km = coe_m.astype(float).copy()
    coe_km[0] /= 1000
    r_j, v_j, info = alignment_at_epoch_to_j2000(
        coe_km, epoch_utc, eop_file, coordinate_epoch=coordinate_epoch_utc)
    rJ2000_m = r_j * 1000
    vJ2000_mps = v_j * 1000
    rvJ2000 = np.vstack((rJ2000_m, vJ2000_mps))
    t = info['OrbitEpochUTC']
    print(f'UTC: {int(t[0]):04d}-{int(t[1]):02d}-{int(t[2]):02d} '
          f'{int(t[3]):02d}:{int(t[4]):02d}:{t[5]:06.3f}')
    print('J2000 r (m)   : [' + '  '.join(f'{x: .6f}' for x in rJ2000_m.ravel()) + ']')
    print('J2000 v (m/s) : [' + '  '.join(f'{x: .6f}' for x in vJ2000_mps.ravel()) + ']')
    print(f"EOP: {info['EOP']['File']}")
    if info['EOP']['IsPredicted']:
        print('该时刻使用预测 EOP 数据。')
    return rJ2000_m, vJ2000_mps, rvJ2000, info


if __name__ == '__main__':
    rJ2000_m, vJ2000_mps, rvJ2000, info = main()
