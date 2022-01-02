import shutil
import glob
import os
import json
import sys
import operator
import math
import copy

import matplotlib
import numpy
import opal_tracking

import xboa.common
import xboa.hit

import config.config_double_triplet_baseline as config
from utils.decoupled_transfer_matrix import DecoupledTransferMatrix
import utils.utilities
import plot

class PlotProbes(object):
    def __init__(self, load_h5, plot_dir, station_in, station_out):
        self.load_h5 = load_h5
        self.plot_dir = plot_dir
        self.station_in = station_in
        self.t_in = [-5, 5]
        self.station_out = station_out
        self.t_out = [110, 120]
        self.var_list = ["r", "r'", "z", "z'"]
        self.weight_list = [1.0, 0.01, 1.0, 0.01]
        self.score_leg = "$\\sqrt{x^2 + y^2 + \\frac{x'^2+y'^2}{0.01^2}}$ [mm]"
        self.bin_width = [0.01, 0.005] # x , y
        self.cmap = "PiYG"
        self.max_z = None


    def get_ke(self, item):
        if item == None:
            return 0.0
        psquared = (item["px"]**2+item["py"]**2+item["pz"]**2)
        mass = xboa.common.pdg_pid_to_mass[2212]
        ke = (psquared + mass**2)**0.5 - mass
        return ke

    def get_plot_list(self):
        data_in = {}
        data_out = {}
        self.plot_list = []
        print("First ", self.load_h5.data[0])
        for item in self.load_h5.data:
            if item["station"] == self.station_in:
                data_in[item["id"]] = item
            elif item["station"] == self.station_out:
                data_out[item["id"]] = item
            else:
                continue
        x_list, y_list, z_list, z2_list, z3_list = [], [], [], [], []

        for id0 in data_in:
            item_in = data_in[id0]
            if id0 in data_out:
                item_out = data_out[id0]
            else:
                item_out = None
            x_list.append(item_in["r"])
            y_list.append(item_in["z"])
            if item_out == None:
                z_list.append(0.1)
            else:
                z_list.append(abs(self.get_ke(item_out)-self.get_ke(item_in)))
                z2_list.append(item_out["z'"]-item_in["z'"])
                z3_list.append(item_out["r'"]-item_in["r'"])
            if id0 < 10:
                print(item_in["r"], item_in["z"], self.get_ke(item_in), self.get_ke(item_out), z_list[-1], item_in["z'"], item_out["z'"])
        return x_list, y_list, z_list, z2_list, z3_list

    def get_bins(self, x_list, bin_width):
        n_bins = int((max(x_list)-min(x_list))/bin_width)+2
        bin_list = [min(x_list)+(i-0.5)*bin_width for i in range(n_bins)]
        return bin_list

    def plot(self):
        x_list, y_list, de_list, zp_list, xp_list = self.get_plot_list()
        x_bin_list = self.get_bins(x_list, self.bin_width[0])
        y_bin_list = self.get_bins(y_list, self.bin_width[1])
        print(x_bin_list)
        print(y_bin_list)
        print(x_list)
        print(y_list)
        print(de_list)
        figure = matplotlib.pyplot.figure()
        axes = figure.add_subplot(1, 1, 1)
        min_z = min(de_list)
        max_z = max(de_list)
        hout = axes.scatter(x_list, y_list,  c=de_list)
        figure.colorbar(hout)
        figure.savefig(os.path.join(self.plot_dir, "dedx.png"))

        figure = matplotlib.pyplot.figure()
        axes = figure.add_subplot(1, 2, 1)
        hout = axes.hist(xp_list)
        axes.text(0.05, 0.95, "RMS r': "+format(numpy.std(xp_list), "4.3g"), transform=axes.transAxes)
        axes.set_xlabel("r'")
        axes = figure.add_subplot(1, 2, 2)
        hout = axes.hist(zp_list)
        axes.text(0.05, 0.95, "RMS z': "+format(numpy.std(xp_list), "4.3g"), transform=axes.transAxes)
        axes.set_xlabel("z'")
        figure.savefig(os.path.join(self.plot_dir, "scattering.png"))

def clear_plot_dir(plot_dir):
    if os.path.exists(plot_dir):
        shutil.rmtree(plot_dir)
    os.makedirs(plot_dir)

def main():
    DecoupledTransferMatrix.det_tolerance = 1.0
    input_dir = "output/arctan_baseline/single_turn_injection/tracking_simulation_v2/foil_test_bumped/grid/"
    plot_dir = input_dir+"/plot_probe_grid/"
    clear_plot_dir(plot_dir)
    plotter = PlotProbes(plot.LoadH5(input_dir+"*FOILPROBE.h5"), plot_dir, 0, 1)
    plotter.max_z = 10
    plotter.plot()

if __name__ == "__main__":
    main()
    matplotlib.pyplot.show(block=False)
    input("Done")
