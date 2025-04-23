from matplotlib import pyplot as plt
from matplotlib import font_manager as fm, rc
import numpy as np

font_legend = fm.FontProperties(family = 'Times New Roman')
tnr = 'Times New Roman'
font = {'family': tnr}

rc('text', usetex=True)
rc('font', **font)

del fm.weight_dict['roman']
# fm._rebuild()

def get_Time (file):
    lines = [list(map(int, line.split())) for line in file.split("\n")[:-1]]
    nv = [line[0] for line in lines]
    time_us = [line[1] for line in lines]
    return nv, time_us
    

def plot_comparison(num_threads, gate = "vanilla"):
    loga = open(f"hyperplonk/logalk {gate} threads {num_threads}.txt").read()
    plk = open(f"hyperplonk/plk {gate} threads {num_threads}.txt").read()
    
    loga_nv, loga_time = get_Time (loga)
    plk_nv, plk_time = get_Time (plk)
    plt.clf()
    plt.title(f"HyperPlonk with Lookup subroutine proof generation - {num_threads} thread(s)", **font, fontweight='bold')
    plt.ylabel("Time (ns)", **font)
    plt.xlabel(r"$\mu$", **font)
    plt.xticks(fontname=tnr)
    plt.yticks(fontname=tnr)
    plt.grid(True)
    plt.plot(loga_nv, loga_time, label="Logalookup")
    plt.plot(plk_nv, plk_time, label="Plookup")
    plt.legend(prop=font_legend)
    
    plt.savefig(f"plot-compare-{num_threads}-threads.png")

def plot_multithread(threads, lk = "logalk", gate = "jellyfish"):
    plt.clf()
    # lk_str = "LogaLookup" if lk == "logalk" else "Plookup"
    # plt.title(f"HyperPlonk with {lk_str} proof generation - Multithread", **font, fontweight='bold')
    plt.ylabel("Time (ns)", **font)
    plt.xlabel(r"$\mu$", **font)
    plt.xticks(fontname=tnr)
    plt.yticks(fontname=tnr)
    plt.grid(True)
    for nt in threads:
        file = open(f"hyperplonk/{lk} {gate} threads {nt}.txt").read()
        nv, time = get_Time (file)
        time = [t * 1000 for t in time]
        plt.plot(nv, time, label=f"{nt} threads")
    plt.legend(prop=font_legend)
    plt.savefig(f"plot-multithread-{lk}.png")

def plot_comparison_lookup_subroutines(num_threads):
    loga = open(f"subroutines/logalk {num_threads} threads.txt").read()
    plk = open(f"subroutines/plk {num_threads} threads.txt").read()
    
    loga_nv, loga_time = get_Time (loga)
    plk_nv, plk_time = get_Time (plk)
    plt.clf()
    # plt.title(f"{num_threads} thread(s)", **font, fontweight='bold')
    plt.ylabel("Time (ns)", **font)
    plt.xlabel(r"$\mu$", **font)
    plt.xticks(fontname=tnr)
    plt.yticks(fontname=tnr)
    plt.grid(True)
    plt.plot(loga_nv, loga_time, label="Logalookup")
    plt.plot(plk_nv, plk_time, label="Plookup")
    plt.legend(prop=font_legend)
    
    plt.savefig(f"plot-compare-lookup-{num_threads}-threads.png")

def plot_multithread_lookup_subroutine(threads, lk = "logalk"):
    plt.clf()
    # lk_str = "LogaLookup" if lk == "logalk" else "Plookup"
    # plt.title(f"{lk_str} subroutine proof generation - Multithread", **font, fontweight='bold')
    plt.ylabel("Time (ns)", **font)
    plt.xlabel(r"$\mu$", **font)
    plt.xticks(fontname=tnr)
    plt.yticks(fontname=tnr)
    plt.grid(True)
    for nt in threads:
        file = open(f"subroutines/{lk} {nt} threads.txt").read()
        nv, time = get_Time (file)
        plt.plot(nv, time, label=f"{nt} threads")
    plt.legend(prop=font_legend)
    plt.savefig(f"plot-multithread-{lk}-subroutine.png")

for nt in [8,16,24,32]:
    plot_comparison(nt)
plot_multithread([8,16,24,32])
for nt in [8,16,24,32]:
   plot_comparison_lookup_subroutines(nt)
plot_multithread_lookup_subroutine([8,16,24,32])