import shutil
import glob
import os
import json
import sys

import numpy
import opal_tracking
import xboa.common
import xboa.hit

import utilities

class PlotProbes(object):
    def __init__(self, file_name_list, plot_dir):
        path = os.path.expandvars("${OPAL_EXE_PATH}/opal")
        mass = xboa.common.pdg_pid_to_mass[2212]
        ref = xboa.hit.Hit.new_from_dict({"pid":2212, "mass":mass, "charge":1.})
        self.tracking = opal_tracking.OpalTracking("", "", ref, file_name_list, path)
        self.plot_dir = plot_dir
        self.probe_data = opal_tracking._opal_tracking.StoreDataInMemory(None)
        

    def load_data(self):
        self.tracking._read_probes(self.probe_data)

    def plot_phase_spaces(self):
        station_list = []
        hit_list = self.probe_data.last[0]
        station_list += [hit['station'] for hit in hit_list]
        station_list = sorted(list(set(station_list)))
        for station in station_list:
            station_hit_list = [hit for hit in hit_list if hit['station'] == station]
            for x_axis , y_axis in [("x", "px"), ("y", "py")]:
                self.plot_phase_space(station_hit_list, x_axis, y_axis)

    def plot_phase_space(self, hit_list, x_axis, y_axis):
        print "\nclosed orbit",
        sys.stdout.flush()
        station = hit_list[0]['station']
        x_list = [hit[x_axis] for hit in hit_list]
        y_list = [hit[y_axis] for hit in hit_list]
        #print x_list, y_list
        data_list = numpy.transpose(numpy.array([x_list, y_list]))
        ell_func = None
        try:
            pass
            #mean, cov = xboa.common.fit_ellipse(data_list, 0.5)
            #ell_func = xboa.common.make_root_ellipse_function(mean, cov)
        except numpy.linalg.linalg.LinAlgError:
            pass
        name = "phase_space_"+str(station).rjust(3, "0")+"_"+x_axis+"_"+y_axis
        canvas = xboa.common.make_root_canvas(name)
        hist, graph = xboa.common.make_root_graph(name, x_list, x_axis, y_list, y_axis)
        time = format(hit_list[0]['t'], "4.4g")
        print station, time
        hist.SetTitle("time "+time+" ns")
        hist.Draw()
        graph.SetMarkerStyle(4)
        graph.Draw("psame")
        if ell_func != None:
            ell_func.Draw("SAME")
        canvas.Update()
        for fmt in ["png"]:
            canvas.Print(self.plot_dir+"/"+name+"."+fmt)

def main():
    a_dir = "output/baseline/tmp/tune/"
    plot_dir = a_dir.split("/tmp/")[0]+"/plot_probe/"
    if os.path.exists(plot_dir):
        shutil.rmtree(plot_dir)
    os.makedirs(plot_dir)
    file_name_list = [a_dir+"RINGPROBE06.loss", a_dir+"RINGPROBE07.loss"]+\
                     glob.glob(a_dir+"CELLPROBE*.loss")
    #file_name_list = [a_dir+"CELLPROBE06.loss"]
    print file_name_list
    plotter = PlotProbes(file_name_list, plot_dir)
    plotter.load_data()
    plotter.plot_phase_spaces()

    #for file_name in :
    #    probe_name = os.path.split(file_name)[1][:-5]
    #    print file_name
    #    data = load_file(file_name)
    #    if len(data) == 0:
    #        continue
    #    plot_phase_space(data, "z", "pz", plot_dir, probe_name)

if __name__ == "__main__":
    main()
    raw_input("Done")
