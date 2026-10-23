"""读取 CelesTrak EOP v1.1 的实测/预测数据，进行逐日线性插值。"""
from pathlib import Path
import numpy as np
from time_utils import valid_epoch, julian_date_midnight


def read_alignment_eop(file, utc):
    utc = valid_epoch(utc)
    file = Path(file).resolve()
    rows, predicted = [], []
    section = ''
    with file.open(encoding='utf-8-sig') as stream:
        for line in stream:
            line = line.strip()
            if line in ('BEGIN OBSERVED', 'BEGIN PREDICTED'):
                section = line.split()[1]
                continue
            if line.startswith('END '):
                section = ''
                continue
            if not section:
                continue
            try:
                values = [float(x) for x in line.split()]
            except ValueError as exc:
                raise ValueError(f'无效 EOP 数据行: {line}') from exc
            if len(values) != 13 or not np.all(np.isfinite(values)):
                raise ValueError(f'无效 EOP 数据行: {line}')
            rows.append(values)
            predicted.append(section == 'PREDICTED')
    rows = np.asarray(rows)
    if rows.size == 0 or np.any(np.diff(rows[:, 3]) <= 0):
        raise ValueError('EOP MJD 必须严格递增且数据不能为空')
    m0 = julian_date_midnight(utc) - 2400000.5
    f = (utc[3] * 3600 + utc[4] * 60 + utc[5]) / 86400
    m = m0 + f
    if not rows[0, 3] <= m <= rows[-1, 3]:
        raise ValueError(f'MJD {m:.9f} 超出 EOP 范围 [{rows[0,3]:.0f}, {rows[-1,3]:.0f}]，不允许外推')
    k = int(np.searchsorted(rows[:, 3], m, side='right') - 1)
    j = k + 1 if m > rows[k, 3] else k
    if rows[k, 3] != m0 or (j != k and rows[j, 3] != m0 + 1):
        raise ValueError('缺少逐日 EOP 记录')
    w = f if j != k else 0.
    v = (1 - w) * rows[k, 4:12] + w * rows[j, 4:12]
    dat = rows[k, 12]
    # DAT 是阶跃量；跨闰秒边界时插值 UT1-TAI。
    du = (1 - w) * (rows[k, 6] - dat) + w * (rows[j, 6] - rows[j, 12]) + dat
    return dict(MJD=m, xp=v[0], yp=v[1], deltaUT1=du, LOD=v[3],
                dPsi=v[4], dEpsilon=v[5], dX=v[6], dY=v[7], deltaAT=dat,
                IsPredicted=bool(predicted[k] or (w > 0 and predicted[j])), File=str(file))
