"""
Plot a single closed orbit (once tracking has finished)
"""
import sys
import copy
import os
import math
import argparse
import h5py
import glob
import shutil


import matplotlib
import ROOT

import xboa.common

import PyOpal.parser
import PyOpal.field
import utils.utilities
import plotting.plot as plot

def plot_azimuthal():
    base_dir = "output/arctan_baseline/single_turn_injection/tracking_simulation_v2/"
    run_dir = os.path.join(base_dir, "foil_test_bumped/grid")
    orbit_folder_list = [(run_dir, "40 mm bumped H$^{+}$"),]
    probe_files = "*PROBE*.h5" #"FOILPROBE_1.h5" #
    track_files = "VerticalSectorFFA-trackOrbi*dat"
    lattice_file = os.path.join(base_dir, "foil_test_bumped/grid/VerticalSectorFFA.tmp")
    angle_domain = [-0.0, 360.0] #[104, 112] #  
    plot_range = [0, 360]
    allowed_events = ["ID1"]

    plot.LoadOrbit.azimuthal_domain = angle_domain
    plot.LoadH5.azimuthal_domain = angle_domain
    test_function = lambda words: words[5] > 1.0 # or math.atan2(float(words[1]), float(words[3])) > 3.0
    orbit_list = []
    for orbit_folder, name in orbit_folder_list:
        h5 = None #plot.LoadH5(os.path.join(orbit_folder, probe_files))
        print(os.path.join(orbit_folder, track_files))
        orbit_file = glob.glob(os.path.join(orbit_folder, track_files))[0]
        print(orbit_folder, track_files, orbit_file)
        orbit = plot.LoadOrbit(orbit_file, allowed_events, test_function)
        orbit_list.append(plot.PlotOrbit(orbit, name, h5))
        print("Loaded", orbit_file)
    plot_fields = plot.PlotFields(plot.GetFields(lattice_file), orbit_list)
    job_name = "closed_orbit_bump"
    figure = plot_fields.azimuthal_fig("Bump", plot_range)
    figure.savefig(os.path.join(base_dir, job_name+".png"))

def plot_beam():
    base_dir = "output/arctan_baseline/isis2_baseline/"
    orbit_folder_list = [(os.path.join(base_dir, "tmp/find_closed_orbits"), "Orbit")]
    probe_files = "*.h5" #"FOILPROBE_1.h5" #
    lattice_file = os.path.join(base_dir, "tmp/find_closed_orbits/VerticalSectorFFA.tmp")
    log_file = os.path.join(base_dir, "tmp/find_closed_orbits/log")
    angle_domain = [-0.0, 360.0] #[104, 112] #  
    plot_range = [108-18, 108+18]
    allowed_events = ["ID1"]
    plot.LoadOrbit.azimuthal_domain = angle_domain
    plot.LoadH5.azimuthal_domain = angle_domain
    test_function = lambda words: False
    orbit_list = []
    log_file = plot.LoadLog(log_file)
    log_file.element_lambda = lambda element: "MAGNET" in element["name"]
    log_file.print()
    for orbit_folder, name in orbit_folder_list:
        orbit_file = os.path.join(orbit_folder, "VerticalSectorFFA-trackOrbit.dat")
        if not os.path.exists(orbit_file):
            orbit_file = os.path.join(orbit_folder, "VerticalSectorFFA-trackOrbit_1.dat")
        orbit = plot.LoadOrbit(orbit_file, allowed_events, test_function)
        orbit_list.append(plot.PlotOrbit(orbit, name, None))
    plot_fields = plot.PlotFields(plot.GetFields(lattice_file), orbit_list)
    plot_fields.polygon_plotter = plot.PlotPolygon(20, 11.2)
    plot_fields.log_plotter = log_file
    plot_fields.n_2d_points = 200
    plot_fields.b0 = 1.51
    job_name = "closed_orbit_bump"
    figure = plot_fields.field_fig("baseline", [0., 0., 0., 0.], [30, 40], [0, 10], False)
    figure.savefig(os.path.join(base_dir, job_name+"_fields2d.png"))
    #figure = matplotlib.pyplot.figure()
    #axes = figure.add_subplot(1, 1, 1)
    #plot_fields.plot_2d(figure, axes, [0.0, 0.0, 0.0, 0.0], [-6, 6], [-6, 6], [-plot_fields.b0, plot_fields.b0], "x", "y", "bz")
    #if plot_fields.polygon_plotter != None:
    #    plot_fields.polygon_plotter.plot(axes)
    #orbit_list[-1].plot_2d(axes, "x", "y", [-5, 5], [-5, 5])
    #axes.set_title("")
    #axes.set_xlabel("x [m]", fontsize=16)
    ##axes.set_ylabel("y [m]", fontsize=16)
    #axes.tick_params(labelsize = 14)
    #matplotlib.pyplot.text(5.10, 5.13, "B$_{z}$ [T]", fontsize=14)
    #figure.savefig(os.path.join(base_dir, job_name+"_fields2d.png"))

def plot_rf():
    base_dir = "output/arctan_baseline/baseline_test_rf_2/"
    orbit_folder_list = [(os.path.join(base_dir, "track_beam_rf_on/grid"), "Orbit")]
    probe_files = "*.h5" #"FOILPROBE_1.h5" #
    lattice_file = os.path.join(base_dir, "track_beam_rf_on/grid/VerticalSectorFFA.tmp")
    log_file = os.path.join(base_dir, "track_beam_rf_on/grid/log")
    angle_domain = [-0.0, 360.0] #[104, 112] #  
    allowed_events = ["ID1"]
    plot.LoadOrbit.azimuthal_domain = angle_domain
    plot.LoadH5.azimuthal_domain = angle_domain
    t0 = 1151.534795/8.0
    r0 = 4.357
    e0 = 0.015
    test_function = lambda words: False
    orbit_list = []
    log_file = plot.LoadLog(log_file)
    log_file.element_lambda = lambda element: "VARIABLE_RF_CAVITY" in element["name"]
    log_file.print()
    for orbit_folder, name in orbit_folder_list:
        orbit_file = os.path.join(orbit_folder, "VerticalSectorFFA-trackOrbit.dat")
        if not os.path.exists(orbit_file):
            orbit_file = os.path.join(orbit_folder, "VerticalSectorFFA-trackOrbit_1.dat")
        orbit = plot.LoadOrbit(orbit_file, allowed_events, test_function)
        orbit_list.append(plot.PlotOrbit(orbit, name, None))
    plot_fields = plot.PlotFields(plot.GetFields(lattice_file), orbit_list)
    plot_fields.polygon_plotter = plot.PlotPolygon(10, 2.8)
    plot_fields.log_plotter = log_file
    plot_fields.n_2d_points = 100
    plot_fields.e0 = e0
    job_name = "baseline"
    figure = plot_fields.rf_fig_2("baseline", [r0*math.cos(math.pi*1.0), r0*math.sin(math.pi*1.0), 0.0, t0], [-4.5, -4.2], [-0.5, 0.5])
    figure.savefig(os.path.join(base_dir, job_name+"_rf.png"))



def main():
    #plot_beam()
    plot_azimuthal()
    #plot_rf()

if __name__ == "__main__":
    main()
    matplotlib.pyplot.show(block = False)
    input("Press <CR> to finish")

