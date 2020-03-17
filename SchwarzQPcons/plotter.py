# Written by Sungho Shin (sungho.shin@wisc.edu), Mihai Anitescu (anitescu@mcs.anl.gov), and Victor Zavala (victor.zavala@wisc.edu)

import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib.ticker import MaxNLocator

plt.rc('text', usetex=True)
plt.rc('font', size=12)
plt.figure(figsize=(6,1.5))
plt.ioff()
plt.rcParams["font.family"] = "serif"


lbl = ["$\\omega=1$","$\\omega=5$","$\\omega=10$"]
clr = [[1,0,0],[0,0,1],[0,1,0]]
mkr = ["o","x","D"]
yl = ["Objective", "Primal Error", "Dual Error"]

for i in range(1):
    tbl=[[0]]*3
    etc=[[0]]*3
    for j in range(3):
        tbl[j] = np.loadtxt("output/tbl-%i-%i.csv" %(i+1,j+1), delimiter=',')
        etc[j] = np.loadtxt("output/etc-%i-%i.csv" %(i+1,j+1), delimiter=',')
    xmax = tbl[0][-1,0]
    for k in range(3):
        for j in range(3):
            plt.plot(tbl[j][:,0],tbl[j][:,k+1],color=clr[j],marker=mkr[j],label=lbl[j],mfc='none')
        if k==0:
            plt.plot([0,xmax],[etc[j][4],etc[j][4]],color=[0,0,0],linestyle='--')
        elif k==1:
            plt.plot([0,xmax],[1e-2,1e-2],color=[0,0,0],linestyle='--')
        elif k==2:
            plt.plot([0,xmax],[1e2,1e2],color=[0,0,0],linestyle='--')
        plt.grid(True,color="lightgrey")
        if k==1:
            plt.legend(loc='lower right',framealpha=1)
        plt.xlabel('Iteration Steps')
        plt.ylabel(yl[k])
        plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
        if k!=0:
            plt.yscale('log')
        plt.grid(True,color="lightgrey")
        plt.xlim(0,60)
        plt.savefig("output/fig-%i-%i.pdf" % (i+1,k+1),bbox_inches="tight")
        plt.clf()
