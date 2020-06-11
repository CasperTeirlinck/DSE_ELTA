import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from collections import OrderedDict
import os

def plotLiftDistribution(y, Cl_range, ClmaxDistr=None, legend=False):
    fig = plt.figure(figsize=(10, 4.5))
    ax1 = fig.add_subplot(111)

    i = 0
    for distr in Cl_range:
        marker = 'o' if i == 0 else 'v' if i == 1 else 's'
        alpha = 2 if i == 0 else 6 if i == 1 else 10
        ax1.plot(y, distr, linewidth=2, color='green', marker=marker, fillstyle='full', markevery=5, label=f'Lift ditribution @ alpha={alpha}')
        if ClmaxDistr: ax1.plot(y, ClmaxDistr(y), linewidth=1.5, linestyle='-.', color='black', label='Local stall limit')
        i += 1

    ax1.axvline(x=0, linewidth=2, color='black')
    ax1.axhline(y=0, linewidth=2, color='black')
    ax1.set_xlabel('Wingspan [m]')
    ax1.set_ylabel('Cl [-]')
    ax1.xaxis.grid(color='black', linestyle='--')
    ax1.yaxis.grid(color='black', linestyle='--')
    if legend: plt.legend(loc='lower center')
    fig.suptitle('Model Lift Distribution', fontsize=16, y=0.97)
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    plt.show()

def readWinglist():
    taper_lst = []
    twist_lst = []
    CLmax_lst = []
    espan_lst = []
    
    with open('Final_Design/aerodynamics/winglist.csv', 'r') as f:
        for line in f.readlines()[1:]:
            columns = line.strip().split(',')
            if len(columns) > 0:
                taper_lst.append(float(columns[0]))
                twist_lst.append(float(columns[1]))
                CLmax_lst.append(float(columns[2]))
                espan_lst.append(float(columns[5]))

    idx = np.argmax(espan_lst)
    bestWing = [taper_lst[idx], twist_lst[idx], CLmax_lst[idx], espan_lst[idx]]

    taper = []
    for i, twist in enumerate(twist_lst):
        if twist == twist_lst[0]:
            taper.append(taper_lst[i])

    CLmax = {}
    for i, twist in enumerate(twist_lst):
        if not twist in CLmax.keys(): CLmax[twist] = []
        CLmax[twist].append(CLmax_lst[i])
    
    espan = {}
    for i, twist in enumerate(twist_lst):
        if not twist in espan.keys(): espan[twist] = []
        espan[twist].append(espan_lst[i])

    return bestWing, taper, OrderedDict(CLmax), OrderedDict(espan)

def plotDesignParams(x, y1, y2, xlabel='', y1label='', y2label=''):
    fig = plt.figure(figsize=(10, 4.5))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    colors = plt.cm.viridis(np.linspace(0, 1, len(y1)))

    for i, twist in enumerate(y1):
        ax1.plot(x, y1[twist], linewidth=2, color=colors[i])
    for i, twist in enumerate(y2):
        ax2.plot(x, y2[twist], linewidth=2, color=colors[i])

    ax1.axvline(x=0, linewidth=2, color='black')
    ax2.axvline(x=0, linewidth=2, color='black')
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(y1label)
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel(y2label)
    ax1.xaxis.grid(color='black', linestyle='--')
    ax1.yaxis.grid(color='black', linestyle='--')
    ax2.xaxis.grid(color='black', linestyle='--')
    ax2.yaxis.grid(color='black', linestyle='--')
    legend = [Line2D([0], [0], color=colors[0], lw=4), Line2D([0], [0], color=colors[int(len(colors)/2)], lw=4), Line2D([0], [0], color=colors[-1], lw=4)]
    plt.legend(legend, [
        f'twist={round(np.degrees(sorted(y1.keys())[0]),0)} deg',
        f'twist={round(np.degrees(sorted(y1.keys())[int(len(y1)/2)]),0)} deg',
        f'twist={round(np.degrees(sorted(y1.keys())[-1]),0)} deg'
    ], loc='lower center')
    fig.suptitle('', fontsize=16, y=0.97)
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    plt.show()