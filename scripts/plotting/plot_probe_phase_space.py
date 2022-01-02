import subprocess
import copy
import shutil
import glob
import os
import json
import sys
import operator
import math

import matplotlib
import numpy
import opal_tracking

import xboa.common
import xboa.hit

import config.config_double_triplet_baseline as config
from utils.decoupled_transfer_matrix import DecoupledTransferMatrix
import utils.utilities

class Cut(object):
    def __init__(self, variable, value, function):
        self.variable = variable
        self.value = value
        self.function = function

    def will_cut(self, hit, decoupled, aa, ref_hit):
        if self.variable == "au" or self.variable == "av":
            index = {"au":1, "av":3}[self.variable]
            return self.function(aa[index], self.value)
        return self.function(hit[self.variable], self.value)

class Cut2(object):
    def __init__(self, variable, value, function, test_hit_list):
        self.rejected = []
        for hit in test_hit_list:
            test_value = hit[variable]
            if function(test_value, value):
                print(hit[-1], test_value, value)
                self.rejected.append(hit[-1])
        #set([hit[-1] for hit in test_hit_list \
        #                                      if function(hit[variable], value)])
        self.rejected = set(self.rejected)
        print("Cutting on", variable, "value", value, "gives:", self.rejected)

    def will_cut(self, hit, decoupled, aa, ref_hit):
        return hit['event_number'] in self.rejected

class RefCut(object):
    def __init__(self, variable, value, function):
        self.variable = variable
        self.value = value
        self.function = function

    def will_cut(self, hit, decoupled, aa, ref_hit):
        return self.function(hit[self.variable]-ref_hit[self.variable], self.value)

class TransmissionCut(object):
    def __init__(self, accept_hit_list):
        self.accepted = set([hit[-1] for hit in accept_hit_list])

    def will_cut(self, hit, decoupled, aa, ref_hit):
        return hit['event_number'] not in self.accepted

class PlotProbes(object):
    def __init__(self, forwards_file_name_list, plot_dir):
        path = os.path.expandvars("${OPAL_EXE_PATH}/opal")
        mass = xboa.common.pdg_pid_to_mass[2212]
        ref = xboa.hit.Hit.new_from_dict({"pid":2212, "mass":mass, "charge":1.})
        file_name_dict = self.setup_stations(forwards_file_name_list) # setup a dict of file_name:station_number
        print(file_name_dict)
        self.tracking = opal_tracking.OpalTracking("", "", ref, file_name_dict, path)
        self.tracking.set_file_format("hdf5")
        self.tracking.name_dict
        self.plot_dir = plot_dir
        my_config = config.Config()
        my_config.tracking["verbose"] = 100
        my_config.tracking["station_dt_tolerance"] = 1000.0
        my_config.tracking["dt_tolerance"] = -1.0
        my_config.tracking["analysis_coordinate_system"] = "azimuthal"
        self.extras = ['t', 'event_number', 'kinetic_energy']
        self.extras_labels = ['time [ns]', 'event number', 'Kinetic Energy [MeV]']
        self.probe_data = opal_tracking._opal_tracking.StoreDataInMemory(my_config)
        self.tracking.pass_through_analysis = self.probe_data
        self.cut_list_1 = []
        self.cut_list_2 = []
        self.tm_dict = {}
        self.ref_dict = {}
        self.stations = {}
        self.max_a4d = 0.040
        self.color_dict = {}
        self.title_text_1 = ""
        self.title_text_2 = ""
        self.fig_list = []
        self.shared_range = []
        self.s = 1
        self.f_size = 20
        self.l_size = 14
        self.ring_tof = 1151.534795
        self.station_by_tof = True
        self.m_index = 1.31
        self.max_amp = []

    def setup_stations(self, file_name_globs):
        file_name_list = []
        for fname in file_name_globs:
            file_name_list += glob.glob(fname)
        endings = [os.path.split(fname)[1] for fname in file_name_list]
        endings = sorted(list(set(endings)))
        station_dict = {}
        ev_counter = {}
        for fname in file_name_list:
            this_ending = os.path.split(fname)[1]
            if this_ending not in ev_counter:
                ev_counter[this_ending] = 0
            ev_counter[this_ending] += 1
            station_dict[fname] = (endings.index(this_ending), ev_counter[this_ending])
        return station_dict


    def load_data(self):
        self.tracking._read_probes()
        for co_params in self.co_param_list:
            self.load_closed_orbits(co_params)
        self.stations = {}
        for hit_list in self.probe_data.last:
            for hit in hit_list:
                if self.station_by_tof:
                    hit["station"] = int(hit["t"]/self.ring_tof)
                if hit["station"] not in self.stations:
                    self.stations[hit["station"]] = 0
                self.stations[hit["station"]] += 1
        print("Loaded following stations:number of hits", json.dumps(self.stations, indent=2))

    def will_cut_1(self, hit, decoupled, aa, ref_hit):
        for cut in self.cut_list_1:
            if cut.will_cut(hit, decoupled, aa, ref_hit):
                return True
        return False

    def will_cut_2(self, hit, decoupled, aa, ref_hit):
        for cut in self.cut_list_2:
            if cut.will_cut(hit, decoupled, aa, ref_hit):
                return True
        return False

    def get_ref_tm(self, hit):
        this_station = hit["station"]
        try:
            ref = self.ref_dict[this_station]
            tm = self.tm_dict[this_station]
        except KeyError:
            ref = self.ref_dict[0]
            tm = self.tm_dict[0]
        ref = copy.deepcopy(ref)
        # dp = p_1.p_0/p_0^2
        dp_over_p = (hit["px"]*ref["px"]+hit["py"]*ref["py"]+hit["pz"]*ref["pz"])/ref["p"]**2-1
        dvert = dp_over_p/self.m_index*1000
        ref["y"] += dvert
        #if hit["station"] < 2:
        #    print("Calculating ref for hit st", hit["station"], "p", hit["p"], "dp_over_p", dp_over_p, "m", self.m_index, "dy", dvert, "y", ref["y"], "hit y", hit["y"])
        return ref, tm


    def get_hits(self, station):
        hit_list_of_lists = self.probe_data.last
        station_hit_list, station_not_list = [], []
        self.cut_list_1.append(Cut("station", station, operator.ne))
        n_cuts = 0
        for hit_list in hit_list_of_lists:
            for hit in hit_list:
                this_station = int(hit["t"]/self.ring_tof) #hit["station"]
                #print("Testing hit with t", hit["t"], "station", this_station, "vs", station)
                if this_station != station:
                    continue
                ref, tm = self.get_ref_tm(hit)
                coupled = [hit[var]-ref[var] for var in ["x", "x'", "y", "y'"]]
                decoupled = tm.decoupled(coupled).tolist()
                action_angle = tm.coupled_to_action_angle(coupled)
                action_angle[1] *= hit["p"]/hit["mass"]
                action_angle[3] *= hit["p"]/hit["mass"]
                coupled = [hit[var] for var in ["x", "x'", "y", "y'"]]
                if self.will_cut_1(hit, decoupled, action_angle, hit_list[0]):
                    n_cuts += 1                   
                    continue
                time = hit["t"]/self.ring_tof
                time = (time - math.floor(time))*self.ring_tof
                energy = hit["kinetic_energy"]
                energy = (hit["p"]-ref["p"])/ref["p"]
                if hit["station"]  < 10:
                    print("  Longitudinal", hit["t"], time, energy)
                if self.will_cut_2(hit, decoupled, action_angle, hit_list[0]):
                    station_not_list.append([this_station]+coupled+decoupled+action_angle+[time, energy]+[hit[e] for e in self.extras]+[ref["y"]])
                else:
                    station_hit_list.append([this_station]+coupled+decoupled+action_angle+[time, energy]+[hit[e] for e in self.extras]+[ref["y"]])
                self.extras_labels += "y$_{ref}$ [mm]"
        if False:
            print("Station", station, "not plotted:", n_cuts, "not accepted:", len(station_not_list), "accepted:", len(station_hit_list))
            for hit in station_hit_list[0:2]:
                print("coupled", end=" ")
                for x in hit[1:5]:
                    print(format(x, "12.6g"), end=" ")
                print("decoupled", end=" ")
                for x in hit[5:9]:
                    print(format(x, "12.6g"), end=" ")
                print("aa", end=" ")
                for x in hit[9:13]:
                    print(format(x, "12.6g"), end=" ")
                print("extras", end=" ")
                for x in hit[13:]:
                    print(format(x, "12.6g"), end=" ")
                print()
        del self.cut_list_1[-1]
        return station_hit_list, station_not_list

    def load_closed_orbits(self, co_params):
        fin = open(co_params["filename"])
        co = json.loads(fin.readline())[0]

        ref_to_bump_station_mapping = co_params["ref_to_bump_station_mapping"]
        for hit in co["ref_track"]:
            #print("hit ", [hit[var] for var in ["station", "x", "px", "y", "py"]])
            #print("seed", [co["seed_hit"][var] for var in ["station", "x", "px", "y", "py"]])
            try:
                bump_station = ref_to_bump_station_mapping[hit["station"]]
            except KeyError: # if hit station not in the list, we ignore it
                continue
            if bump_station in self.ref_dict:
                continue
            mass = xboa.common.pdg_pid_to_mass[2212]
            hit.update({"pid":2212, "mass":mass, "charge":1.})
            hit = xboa.hit.Hit.new_from_dict(hit)
            self.ref_dict[bump_station] = hit
            self.ref_dict[bump_station+100] = hit
        tm = [row[1:5] for row in co["tm"]]
        tm = DecoupledTransferMatrix(tm, True)
        for var in ref_to_bump_station_mapping.values():
            self.tm_dict[var] = tm
        #tm.print_tests()

    def plot_phase_spaces(self):
        station_list = sorted(list(self.stations.keys()))
        station_list = station_list[0:100:1]#+station_list[60::10]#+station_list[0:11:1]+station_list[60::10]+station_list[-2:-1:1]
        station_list = set(station_list)
        for station in station_list:
            h_list, n_list = self.get_hits(station)
            if len(self.color_dict) == 0:
                self.build_colors(h_list)
                self.s = 2*utils.utilities.matplot_marker_size(h_list)
            if len(h_list) == 0:
                print("No hits for station", station)
                continue
            z_axis = None
            xlim, ylim = None, None
            self.max_amp.append(max([h[12] for h in h_list]))
            for name, hit_list, not_list in [("", h_list, n_list)]:#, ("not-", [], n_list), ("hit-", h_list, [])]:
                figure, axis_list = utils.utilities.setup_da_figure(False)
                self.plot_phase_space(axis_list[0], hit_list, not_list, 1, 3, z_axis, station)
                self.plot_phase_space(axis_list[1], hit_list, not_list, 1, 2, z_axis, station)
                self.plot_phase_space(axis_list[2], hit_list, not_list, 3, 4, z_axis, station)
                self.plot_phase_space(axis_list[3], hit_list, not_list, 5, 7, z_axis, station)
                self.plot_phase_space(axis_list[4], hit_list, not_list, 5, 6, z_axis, station)
                self.plot_phase_space(axis_list[5], hit_list, not_list, 7, 8, z_axis, station)
                self.plot_phase_space(axis_list[6], hit_list, not_list, 10, 12, z_axis, station)
                self.plot_phase_space(axis_list[7], hit_list, not_list, 9, 10, z_axis, station)
                self.plot_phase_space(axis_list[8], hit_list, not_list, 11, 12, z_axis, station)
                if name == "":
                    self.title_text_2 = "showing good and bad hits"
                elif name == "not-":
                    self.title_text_2 = "showing only bad hits"
                elif name == "hit-":
                    self.title_text_2 = "showing only good hits"
                else:
                    self.title_text_2 = name
                #self.suptitle(figure, station, hit_list, not_list)
                station_str = str(station).rjust(3, '0')
                figure.suptitle("Station: "+str(station))#+"\nwith "+str(len(hit_list))+" good hits")
                figure.savefig(self.plot_dir+"/transverse_station-"+name+station_str+".png")
                self.fig_list.append(figure)
                matplotlib.pyplot.close(figure)
                figure = matplotlib.pyplot.figure(figsize=(10, 6))
                figure.suptitle("Station: "+str(station))
                axes = figure.add_subplot(1, 1, 1,  position=[0.15, 0.2, 0.3, 0.7])
                self.plot_phase_space(axes, hit_list, not_list, 13, 14, z_axis, station)
                axes = figure.add_subplot(1, 1, 1,  position=[0.6, 0.2, 0.3, 0.7])
                self.plot_phase_space(axes, hit_list, not_list, 3, 14, z_axis, station)
                axes.plot(axes.get_ylim(), axes.get_ylim())
                figure.savefig(self.plot_dir+"/longitudinal_station-"+name+station_str+".png")
                self.fig_list.append(figure)
                matplotlib.pyplot.close(figure)


    def get_lim(self, var):
        lim_dict = {
            1:[4300, 4400],
            2:[-0.01, 0.01],
            3:[-70, 70],
            4:[-0.01, 0.01],
            5:[-70, 70],
            6:[-0.01, 0.01],
            7:[-100, 100],
            8:[-0.01, 0.01],
            9:[-math.pi, math.pi],
            10:[-0.0001, 0.04],
            11:[-math.pi, math.pi],
            12:[-0.0001, 0.04],
            13:[0.0, self.ring_tof],
            14:[-0.025, 0.025],
        }
        lim_dict[18] = lim_dict[3]
        if var in lim_dict:
            return lim_dict[var]
        return None

    def plot_phase_space(self, axes, hit_list, not_list, x_axis, y_axis, z_axis, station):
        x_list = [hit[x_axis] for hit in hit_list]
        y_list = [hit[y_axis] for hit in hit_list]
        if z_axis != None:
            z_list = [hit[z_axis] for hit in hit_list]
        else:
            z_list = [self.color(hit) for hit in hit_list]

        x_not_list = [hit[x_axis] for hit in not_list]
        y_not_list = [hit[y_axis] for hit in not_list]
        labels = {
            "kinetic_energy":"KE [MeV]",
            "x":"x [mm]",
            "y":"y [mm]",
            "px":"p$_x$ [MeV/c]",
            "py":"p$_y$ [MeV/c]",
            "x'":"x'",
            "y'":"y'",
            "t":"time [ns]",
            0:"station",
            1:"x [mm]",
            2:"p$_x$ [MeV/c]",
            3:"y [mm]",
            4:"p$_y$ [MeV/c]",
            5:"u",
            6:"p$_u$",
            7:"v",
            8:"p$_v$",
            9:"$\\phi_u$",
            10:"Norm. A$_u$ [mm]",
            11:"$\\phi_v$",
            12:"Norm. A$_v$ [mm]",
            13:"$\\delta$t [ns]",
            14:"dp/p",
            18:"y$_{ref}$ [mm]"
        }
        for i, label in enumerate(self.extras_labels):
            labels[i+15] = label
        if True:
            labels[2] = "x'"
            labels[4] = "y'"
            labels[6] = "u'"
            labels[8] = "v'"
        name = "phase_space_"+str(station).rjust(3, "0")+"_"+str(x_axis)+"_"+str(y_axis)
        if len(not_list):
            scat = axes.scatter(x_not_list, y_not_list, c='limegreen', edgecolors=None, s=self.s)
        if len(hit_list):
            scat = axes.scatter(x_list, y_list, c=z_list, edgecolors=None, s=self.s)
        axes.set_xlabel(labels[x_axis], fontsize=self.f_size)
        axes.set_ylabel(labels[y_axis], fontsize=self.f_size)
        axes.grid(True)
        x_lim, y_lim = self.get_lim(x_axis), self.get_lim(y_axis)
        if x_lim:
            axes.set_xlim(x_lim)
        if y_lim:
            axes.set_ylim(y_lim)
        axes.tick_params(labelsize = self.l_size)
        if z_axis != None:
            axes.get_figure().colorbar(scat)

    def color(self, hit):
        return 0 #self.color_dict[hit[16]] # event number

    def color_hit(self, hit):
        au = hit[10]
        av = hit[12]
        hue = math.atan2(au, av)/(math.pi/2)*(1./3.)+2./3. # range between 0 and 2/3 i.e. green through red through blue
        saturation = min(1, 0.2+0.8*(au+av)/self.max_a4d)
        value = 0.9
        rgb = matplotlib.colors.hsv_to_rgb([hue, saturation, value]) 
        return rgb

    def build_colors(self, hit_list):
        self.color_dict = {}
        a4d_list = [hit[10]+hit[12] for hit in hit_list]
        self.max_a4d = max(a4d_list)
        for hit in hit_list:
            self.color_dict[hit[16]] = self.color_hit(hit)

    def suptitle(self, figure, station, hit_list, not_hit_list):
        title = "Station "+str(station)+" Good hits: "+str(len(hit_list))+"\n"
        t_list = [hit[-2] for hit in hit_list]
        title += " $\\bar{t}$: "+format(numpy.mean(t_list), "6.2f")+" ns"
        title += " $\\sigma(t)$: "+format(numpy.std(t_list), "6.2f")+" ns"
        if self.title_text_1:
            title += "\n"+self.title_text_1+" "+self.title_text_2

        figure.suptitle(title)


    def movie(self):
        here = os.getcwd()
        os.chdir(self.plot_dir)
        #mencoder mf://turn*.png -mf w=800:h=600:fps=5:type=png -ovc lavc -lavcopts vcodec=msmpeg4:mbd=2:trell -oac copy -o injection.avi
        try:
            output = subprocess.check_output(["mencoder",
                                    "mf://transverse_*.png",
                                    "-mf", "w=800:h=600:fps=5:type=png",
                                    "-ovc", "lavc",
                                    "-lavcopts", "vcodec=msmpeg4:vbitrate=2000:mbd=2:trell",
                                    "-oac", "copy",
                                    "-o", "transverse_injection.avi"])
        except:
            print("Transverse movie failed")
        try:
            output = subprocess.check_output(["mencoder",
                                    "mf://longitudinal_*.png",
                                    "-mf", "w=800:h=600:fps=5:type=png",
                                    "-ovc", "lavc",
                                    "-lavcopts", "vcodec=msmpeg4:vbitrate=2000:mbd=2:trell",
                                    "-oac", "copy",
                                    "-o", "longitudinal_injection.avi"])
        except:
            print("Longitudinal movie failed")
        os.chdir(here)


def main():
    DecoupledTransferMatrix.det_tolerance = 1.0
    glob_dir = "output/arctan_baseline/single_turn_injection/tracking_simulation_ref-only-hor/26-mm/rf-and-foil"
    delta = []
    max_amp = []
    for a_dir in glob.glob(glob_dir):
        #a_dir = "output/arctan_baseline/baseline_test_rf_2/track_beam_rf_on"
        plot_dir = a_dir+"/plot_probe_phase_foil/"
        if os.path.exists(plot_dir):
            shutil.rmtree(plot_dir)
        os.makedirs(plot_dir)
        for probe in ["4"]:
            for_glob_name = a_dir+"/grid/RINGPROBE0"+str(probe)+".h5"
            #for_glob_name = a_dir+"/grid/FOILPROBE_"+str(probe)+".h5"
            forwards_file_name_list = glob.glob(for_glob_name)
            plotter = PlotProbes(forwards_file_name_list, plot_dir)
            plotter.co_param_list = [{
                "filename":os.path.join("output/arctan_baseline/baseline/closed_orbits_cache"),
                "ref_to_bump_station_mapping":dict([(i,i) for i in range(1001)]),
            },]
            try:
                plotter.load_data()
            except IOError:
                print("IOError trying", for_glob_name)
                raise
            #plotter.stations = set([0])
            cut_station = 25
            cut_hits = plotter.get_hits(cut_station)[0]
            # events that fail cut_list_1 are not plotted at all
            #plotter.cut_list_1.append(Cut("pz", 0., operator.lt))
            #plotter.cut_list_1.append(Cut("au", 0.5, operator.gt))
            #plotter.cut_list_1.append(Cut("av", 0.5, operator.gt))
            # events that fail cut_list_2 are put in the "not" list and plotted as hollow dots
            #plotter.cut_list_2.append(TransmissionCut(cut_hits))
            #plotter.cut_list_2.append(Cut2(10, 0.20, operator.gt, cut_hits))
            #plotter.cut_list_2.append(Cut2(12, 0.20, operator.gt, cut_hits))
            #plotter.title_text_1 = "Cut at turn "+str(cut_station)
            plotter.plot_phase_spaces()
            plotter.movie()
            delta.append(float(a_dir.split("ref-only-")[1].split("-mm")[0]))
            max_amp.append(max(plotter.max_amp[-9:]))
    figure = matplotlib.pyplot.figure()
    axes = figure.add_subplot(1, 1, 1)
    axes.semilogy(delta, max_amp)
    axes.set_xlabel("nominal bump [mm]")
    axes.set_ylabel("Norm. A$_v$ [mm]")
    figure.savefig(a_dir+"/../../amplitude_vs_bump.png")

if __name__ == "__main__":
    main()
    matplotlib.pyplot.show(block=False)
    input("Done")
