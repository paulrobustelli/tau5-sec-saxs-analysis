"""Reproduce DENSS reconstruction and alignment; run from this directory.
Install requirements.txt first. 20 independently initialized maps per sample.
DENSS records each automatically generated random seed in its run log.
"""
import os,sys,subprocess
from pathlib import Path
os.chdir(Path(__file__).resolve().parent)
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
for sample,dmax in [('WT',91.7257817),('AA',166.7131223)]:
 subprocess.run([sys.executable,'-m','denss.scripts.denss_all','-f',sample+'.dat','-d',str(dmax),'-m','SLOW','-nm','20','-j','4','--plot_off','--write_freq','10000','-o',sample+'_denss'],check=True)

for sample,dmax in [('WT',91.7257817),('AA',166.7131223)]:
 subprocess.run([sys.executable,'-m','denss.scripts.denss_refine','-f',sample+'.dat','-d',str(dmax),'-rho',sample+'_denss/'+sample+'_denss_avg.mrc','--seed','20260923','--plot_off','--write_freq','10000','-o',sample+'_refined'],check=True)
