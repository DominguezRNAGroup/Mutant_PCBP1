from pathlib import Path
import mdtraj as md
import pickle
from math import log10, floor

clusterp = Path('rmsd-clustering-2025-12-02')
for cp in clusterp.iterdir():
    print(cp)
    if cp.is_dir():
        centersp = cp/'centers.pickle'
        cendir = cp/'centers'
        cendir.mkdir(exist_ok=True)
        centers = pickle.loads(centersp.read_bytes())
        centers[0].save(str(cendir/'center0.pdb'))
        ncens = len(centers)
        pad = int(floor(log10(ncens))) + 1
        print(cp.name, ncens)
        merged_cens = md.join(centers)
        merged_cens.superpose(merged_cens)
        merged_cens.save(str(cendir/'centers.dcd'))
        # for i, cen in enumerate(centers):
            # cen.save(cendir/f'center-{i:0>{pad}}.pdb')
