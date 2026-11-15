import shutil
from common import *
pairs=[(f'{HOME}/jobs/ncasi_p2_20260908/figs_notitle/fia-vs-treemap-divergence_DATA_2026-09-07.csv','fig3_treemap_phase2.csv','Verbatim copy (already clean, aggregate model-pair contrast).'),
       (f'{HOME}/perseus_me/tracks/multistate/tracks/multistate/ME_basis_crosswalk.csv','ME_basis_crosswalk.csv','Verbatim copy of the multistate release crosswalk (multistate vs Maine a3 basis).')]
for src,dst,note in pairs:
    shutil.copyfile(src,f'{DATA}/{dst}'); assert md5(src)==md5(f'{DATA}/{dst}')
    record(dst,[src],__file__,note)
    status(f'step3/6 copied {dst}')
