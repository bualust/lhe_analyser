import awkward as ak
import pylhe
import math
import boost_histogram as bh
import matplotlib.pyplot as plt
import mplhep as hep
from matplotlib.gridspec import GridSpec
import os
import argparse
import yaml

def main():
    parser = argparse.ArgumentParser(description="Process some integers.")
    parser.add_argument(
        "--n",
        dest="nevents",
        type=int,
        required=False,
        metavar="1000",
        default=-1,
        help="Number of events to plot",
    )
    parser.add_argument(
        "--o",
        dest="outdir",
        type=str,
        required=True,
        metavar="output_plots",
        help="Name of dir containing plots",
    )
    parser.add_argument(
        "--i",
        dest="indir",
        type=str,
        required=False,
        metavar="lhe-files/",
        help="Name of dir containing lhe files",
        default="lhe-files/",
    )
    parser.add_argument(
        "--c",
        dest="config",
        type=str,
        required=True,
        metavar="config.yaml",
        help="Configuration file with processes to compare",
    )

    args = parser.parse_args()

    allgood(f"Running of {args.nevents} and saving plots in {args.outdir}")

    with open(args.config) as yaml_file:
        yaml_config = yaml.load(yaml_file, Loader=yaml.FullLoader)
    processes = yaml_config

    process_plot = {}
    variables = {"M(bb)":[100,95,150],
                 "M(ee)":[100,60,120],
                 "pT(b)":[100,0,200],
                 "pT(anti-b)":[100,0,200],
                 "DeltaR(bb)":[100,0,4],
                 "DeltaR(bb,ee)":[100,2.5,6],
                 "DeltaR(b,ee)":[100,2.5,6],
                 "DeltaR(anti-b,ee)":[100,2.5,6],
                 "pT(b) x spin":[100,-200,200],
                 "pT(anti-b) x spin":[100,-200,200],
                 "spin(b)":[6,-3,3],
                 "spin(anti-b)":[6,-3,3],
                 "DeltaAngle(bb)":[100,0,4],
                 "DeltaAngle(bxb,ee)":[100,0,4],
                 "weight":[100,0,1]
                 }

    for key, value in processes.items():
        lhe_file = f"{args.indir}{value}unweighted_events.lhe"
        allgood(f"Reading file:\n{lhe_file}")
        allgood(f"Number of events: {pylhe.read_num_events(lhe_file)}")

        pylhe_file = pylhe.read_lhe_with_attributes(lhe_file)
        events = pylhe.to_awkward(pylhe_file)

        plots = {}
        if args.nevents != -1:
            particles = events.particles[:args.nevents]
            plots["weight"] = events.eventinfo.weight[:args.nevents]
        else:
            particles = events.particles
            plots["weight"] = events.eventinfo.weight

        stable_parts = stable_particles(particles)

        plots["M(bb)"] = mass_bb(stable_parts)
        plots["M(ee)"] = mass_ee(stable_parts)
        plots["pT(b)"] = pt_b(stable_parts,1)
        plots["pT(anti-b)"] = pt_b(stable_parts,-1)
        plots["DeltaR(bb)"] = deltar_bb(stable_parts)
        plots["DeltaR(bb,ee)"] = deltar_bb_ee(stable_parts,0)
        plots["DeltaR(b,ee)"] = deltar_bb_ee(stable_parts,1)
        plots["DeltaR(anti-b,ee)"] = deltar_bb_ee(stable_parts,-1)
        plots["pT(b) x spin"] = pt_b_times_spin(stable_parts,1)
        plots["pT(anti-b) x spin"] = pt_b_times_spin(stable_parts,-1)
        plots["spin(b)"] = b_spin(stable_parts,1)
        plots["spin(anti-b)"] = b_spin(stable_parts,-1)
        plots["DeltaAngle(bb)"] = deltaangle_bb(stable_parts)
        plots["DeltaAngle(bxb,ee)"] = cross_bb_ee(stable_parts)

        process_plot[key] = plots

    for var,binning in variables.items():
        if var=="weight":
            continue
        fig = plt.figure(figsize=(7.0,7.0),dpi=400)
        gs = GridSpec(2,1, height_ratios=[6,1],hspace=0.1)
        main = fig.add_subplot(gs[0])
        ratio= fig.add_subplot(gs[1],sharex=main)

        info(f"{var} had binning {binning}")
        for proc, plots_dic in process_plot.items():
            h_ax = bh.axis.Regular(binning[0],binning[1],binning[2])
            hist = bh.Histogram(h_ax,storage=bh.storage.Weight())
            if args.nevents!=-1:
                hist.fill(plots_dic[var], weight=plots_dic["weight"][:args.nevents])
            else:
                hist.fill(plots_dic[var], weight=plots_dic["weight"])
            hep.histplot(hist,yerr=False,histtype='step',
                     ax=main,density=True,label=proc)

        plt.xlabel(var)
        ymin, ymax = main.axes.get_ylim()
        main.set_ylim(ymin,ymax+(ymax-ymin)*0.3)
        ratio.set_ylim(0.0,2.0)
        main.axes.xaxis.set_visible(False)
        main.legend(loc='upper right')
        os.system('mkdir -p '+args.outdir)
        ratio.axhline(y=1.0, color='grey',linestyle='--')
        fig.savefig(f"{args.outdir}/{var}.pdf")

def stable_particles(particles):
    return particles[particles.status == 1]

def mass_bb(stable_parts):
    b_quark = stable_parts[stable_parts.id == 5]
    antib_quark = stable_parts[stable_parts.id == -5]
    mbb = (b_quark.vector[:]+antib_quark.vector[:]).mass
    return ak.flatten(mbb)

def deltar_bb_ee(stable_parts,which):
    b_quark = stable_parts[stable_parts.id == 5]
    antib_quark = stable_parts[stable_parts.id == -5]
    el = stable_parts[stable_parts.id == 11]
    antiel = stable_parts[stable_parts.id == -11]
    bb_pair = b_quark.vector[:].add(antib_quark.vector[:])
    ee_pair = el.vector[:].add(antiel.vector[:])
    deltar = bb_pair.deltaR(ee_pair)
    deltarb = b_quark.vector[:].deltaR(ee_pair)
    deltarantib = antib_quark.vector[:].deltaR(ee_pair)
    if which == 0:
        return ak.flatten(deltar)
    elif which == 1:
        return ak.flatten(deltarb)
    elif which == -1:
        return ak.flatten(deltarantib)

def deltar_bb(stable_parts):
    b_quark = stable_parts[stable_parts.id == 5]
    antib_quark = stable_parts[stable_parts.id == -5]
    deltar = (b_quark.vector[:].deltaR(antib_quark.vector[:]))
    return ak.flatten(deltar)

def cross_bb_ee(stable_parts):
    b_quark = stable_parts[stable_parts.id == 5]
    antib_quark = stable_parts[stable_parts.id == -5]
    b_quark_3d = b_quark.vector[:].to_Vector3D()
    antib_quark_3d = antib_quark.vector[:].to_Vector3D()
    cross_bb = (b_quark_3d.cross(antib_quark_3d))

    el = stable_parts[stable_parts.id == 11]
    antiel = stable_parts[stable_parts.id == -11]
    el_3d = el.vector[:].to_Vector3D()
    antiel_3d = antiel.vector[:].to_Vector3D()
    ee_pair = el_3d.add(antiel_3d)

    deltaangle = ee_pair.deltaangle(cross_bb)
    return ak.flatten(deltaangle)

def deltaangle_bb(stable_parts):
    b_quark = stable_parts[stable_parts.id == 5]
    antib_quark = stable_parts[stable_parts.id == -5]
    deltaangle = (b_quark.vector[:].deltaangle(antib_quark.vector[:]))
    return ak.flatten(deltaangle)

def mass_ee(stable_parts):
    el = stable_parts[stable_parts.id == 11]
    antiel = stable_parts[stable_parts.id == -11]
    mee = (el.vector[:]+antiel.vector[:]).mass
    return ak.flatten(mee)

def pt_b(stable_parts,charge):
    b_quark = stable_parts[stable_parts.id == 5*charge]
    pt_b = (b_quark.vector[:]).pt
    return ak.flatten(pt_b)

def b_spin(stable_parts,charge):
    b_quark = stable_parts[stable_parts.id == 5*charge]
    spin_b = b_quark.spin[:]
    return ak.flatten(spin_b)

def pt_b_times_spin(stable_parts,charge):
    b_quark = stable_parts[stable_parts.id == 5*charge]
    pt_b = ((b_quark.vector[:]).pt)*b_quark.spin[:]
    return ak.flatten(pt_b)

class bcolors:
    INFO = '\033[95m' #PURPLE
    OK = '\033[92m' #GREEN
    WARNING = '\033[93m' #YELLOW
    FAIL = '\033[91m' #RED
    RESET = '\033[0m' #RESET COLOR

def warn(message): print(bcolors.WARNING+"=== "+message+" ==="+bcolors.RESET)
def fail(message): print(bcolors.FAIL+"=== "+message+" ==="+bcolors.RESET)
def allgood(message): print(bcolors.OK+"=== "+message+" ==="+bcolors.RESET)
def info(message): print(bcolors.INFO+"=== "+message+" ==="+bcolors.RESET)

if __name__ == "__main__":
    main()
