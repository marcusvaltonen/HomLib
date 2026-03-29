import glob

import matplotlib.pyplot as plt
import numpy as np

solvers_one_sided = {
    "lomsac_fitzgibbon_cvpr_2001_one_sided": 'fitzgibbon_one_sided',
    "lomsac_nakano_icpr_2025_one_sided": 'nakano_one_sided',
    "lomsac_wadenback_2025_one_sided": 'wadenback_one_sided',
}

solvers_two_sided_equal = {
    "lomsac_kukelova_cvpr_2015_two_sided_equal": 'kukelova_two_sided_equal',
    "lomsac_kukelova_cvpr_2015_two_sided_equal_6pt": 'kukelova_two_sided_6pt_equal',
    "lomsac_fitzgibbon_cvpr_2001_two_sided_equal": 'fitzgibbon_equal',
    "lomsac_wadenback_2025_two_sided_equal": 'wadenback_two_sided_equal',
}


solvers_two_sided = {
    "lomsac_kukelova_cvpr_2015_two_sided": 'kukelova_two_sided',
    "lomsac_kukelova_cvpr_2015_two_sided_6pt": 'kukelova_two_sided_6pt',
    "lomsac_wadenback_2025_two_sided": 'wadenback_two_sided',
}


solvers = {
    #**solvers_one_sided,
    #**solvers_two_sided_equal,
    **solvers_two_sided,
}


container = {}

for s, ss in solvers.items():
    files = glob.glob(f'hpatches_time_mAA_{s}_*.npy')
    if not files:
        print('NO FILES!!!')
    data = []
    for f in files:
        if "equal" in f:
            continue
        if (s == "lomsac_kukelova_cvpr_2015_two_sided_equal" or s == "lomsac_kukelova_cvpr_2015_two_sided") and '6pt' in f:
            continue
        d = np.loadtxt(f)
        data.append(d)

    data = np.array(data)
    i = data[:, 0].argsort(axis=0)
    data[:,0] *= 1e-6  # Convert to ms
    data = data[i, :]
    np.savetxt(f'hpatches_time_mAA_{ss}.csv', data, delimiter=',')

