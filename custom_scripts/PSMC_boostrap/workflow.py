from workflow_targets import *
from gwf import *
gwf = Workflow()

ROOT_PATH="/home/jilong/EcoGenetics/people/jilong/psmc"
for sp in ["DUM","TEN","SAR","BIC","MIM","LIN"]:
    for chrom_type in ["X","A"]:
        if chrom_type == "A":
            mu = 5e-09
        elif chrom_type == "X":
            mu = 3.8e-09
        g = 1
        n = 100
        run_n_bootstrap_psmc(gwf,ROOT_PATH,sp,chrom_type,n,mu,g)
