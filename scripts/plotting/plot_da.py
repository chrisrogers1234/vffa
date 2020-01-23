import copy
import math
import glob
import os
import json
import sys

import ROOT
import numpy
import scipy.spatial
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot

import xboa.common as common
from utils import utilities
from xboa.hit import Hit
from xboa.algorithms.tune import DPhiTuneFinder
from utils.decoupled_transfer_matrix import DecoupledTransferMatrix

ROOT_OBJECTS = []



class DAPlotter(object):
    def __init__(self, file_name):
        self.load_data(file_name)
        self.var_list = ['x', 'px', 'y', 'py']
        self.plot_dir = "./"
        self.formats = ["png"]
        DecoupledTransferMatrix.det_tolerance = 1.

    def load_data(self, file_name):
        fin = open(file_name)
        self.data = []
        for line in fin.readlines():
            try:
                self.data.append(json.loads(line))
            except:
                continue
        print("Loaded", len(self.data), "lines")

    def polar(self, point):
        polar = [math.atan2(point[1], point[0]),
                (point[0]**2+point[1]**2)**0.5]
        return polar


    def centre(self, points_cartesian):
        mean = [
            numpy.mean([point[0] for point in points_cartesian]),
            numpy.mean([point[1] for point in points_cartesian])
        ]
        points_cartesian -= mean
        return points_cartesian

    def calculate_da(self, row_data, da_key, n_points):
        return 0.
        axis1, axis2 = {"x_da":("x", "x'"), "y_da":("y", "y'")}[da_key]
        points = []
        mass = common.pdg_pid_to_mass[2212]
        hit_list = []
        for i, tmp_list in enumerate(row_data[da_key]):
            print(len(tmp_list), tmp_list[0])
            # hit_list[0] is the x/y offset in mm; hit_list[1] is the corresponding
            # set of probe hits
            if len(tmp_list[1]) < n_points:
                break
            hit_list = tmp_list
        hit_list = [Hit.new_from_dict(hit_dict) for hit_dict in hit_list[1]]
        points = numpy.array([[hit[axis1], hit[axis2]] for hit in hit_list])
        x_list = [hit[axis1] for hit in hit_list]
        print(min(x_list), max(x_list))
        y_list = [hit[axis2] for hit in hit_list]
        print(min(y_list), max(y_list))
        pz = hit_list[0]['pz']
        print(math.pi*(max(x_list)-min(x_list))*(max(y_list)-min(y_list))/4.*mass/pz)
        if len(points):
            geometric_acceptance = get_area(points)/math.pi
            # beta gamma * geometric emittance
            normalised_acceptance = geometric_acceptance*pz/mass
            return normalised_acceptance
        else:
            return 0.

    def get_area(self, points_cartesian):
        points_cartesian = self.centre(points_cartesian)
        points_polar = [self.polar(point) for point in points_cartesian]
        points_polar = numpy.array(sorted(points_polar))
        area = 0.
        for i, point_1 in enumerate(points_polar[1:]):
            point_0 = points_polar[i]
            # triangle with base length r0 and height r1 sin phi
            phi = point_1[0] - point_0[0]
            h = point_1[1]*math.sin(phi)
            delta_area = 0.5*point_0[1]*h
            area += delta_area
        phi_first = 2.*math.pi+points_polar[0, 0]-points_polar[-1, 0]
        delta_area = 0.5*points_polar[0, 1]*math.sin(phi_first)*points_polar[-1, 1]
        area += delta_area
        return area

    def plot_da(self, max_n_points, plot_dir, acceptance):
        variables = utilities.get_substitutions_axis(self.data, "substitutions")
        plot_dir = os.path.join(plot_dir, "plot_da")
        utilities.clear_dir(plot_dir)
        self.plot_dir = plot_dir
        for index, item in enumerate(self.data):
            for key in variables:
                variables[key] = item['substitutions'][key]
            for da_key in 'y_da', 'x_da': 
                if da_key not in list(item.keys()):
                    print("Did not find", da_key, "in keys", list(item.keys()), "... skipping")
                    continue
                for axis1, axis2 in ('x', "x'"), ('y', "y'"):
                    canvas = self.plot_one_da(item, da_key, axis1, axis2, max_n_points, variables, acceptance)
                    plot_name = "da_"+str(index)+"_"+da_key+"_"+axis1+"_"+axis2
                    for format in ["eps", "png", "root"]:
                        canvas.Print(plot_dir+"/"+plot_name+"."+format)

                    continue
                    if axis1[0] not in da_key:
                        continue
                    canvas, chol_canvas = plot_one_phi(item, da_key, axis1, axis2, max_n_points, variables)
                    plot_name = "phi_"+str(index)+"_"+da_key+"_"+axis1+"_"+axis2
                    for format in ["eps", "png", "root"]:
                        canvas.Print(plot_dir+"/"+plot_name+"."+format)
                        chol_canvas.Print(plot_dir+"/"+plot_name+"_cholesky."+format)


                self.plot_decoupled_da(item, da_key, axis1, axis2, max_n_points, variables, acceptance)

    def make_global_plot(self, min_n_points, plot_dir):
        axis_candidates = utilities.get_substitutions_axis(self.data, "substitutions")
        my_axes_x = dict([(key, []) for key in axis_candidates])
        da = {'y_da':[], 'x_da':[]}
        for item in self.data:
            for da_key in da:
                if da_key not in list(item.keys()):
                    print("Did not find", da_key, "in keys", list(item.keys()), "... skipping")
                    continue
                da[da_key].append(self.calculate_da(item, da_key, min_n_points)*1e3)
                print(da_key, da[da_key][-1])
            for ax_key in my_axes_x:
                my_axes_x[ax_key].append(item['substitutions'][ax_key])
        for ax_key in my_axes_x:
            ordinate = my_axes_x[ax_key]
            x_name = utilities.sub_to_name(ax_key)
            canvas_name = x_name+' vs da'
            canvas = common.make_root_canvas(canvas_name)
            canvas.SetLogy()
            x_range = [min(ordinate), max(ordinate)]
            y_range = [min(da['x_da']+da['y_da']), max(da['x_da']+da['y_da'])]
            y_range[0] = max(y_range[0]/2, 0.01)
            hist, graph = common.make_root_graph('axes', x_range, x_name, y_range, "DA [#pi #mum]", ymin=y_range[0], ymax=y_range[1]*3)
            hist.Draw()
            if len(da['x_da']) == len(ordinate):
                hist_x, graph_x = common.make_root_graph('horizontal da', ordinate, "", da['x_da'], "")
                graph_x.SetMarkerStyle(24)
                graph_x.SetMarkerColor(2)
                graph_x.Draw("SAME P")
            if len(da['y_da']) == len(ordinate):
                hist_y, graph_y = common.make_root_graph('vertical da', ordinate, "", da['y_da'], "")
                graph_y.SetMarkerStyle(26)
                graph_y.SetMarkerColor(4)
                graph_y.Draw("SAME P")
            #legend = common.make_root_legend(canvas, [graph_x, graph_y])
            #legend.Draw()
            canvas.Update()
            canvas_name = canvas_name.replace(" ", "_")
            for format in "eps", "png", "root":
                canvas.Print(plot_dir+"/"+canvas_name+"."+format)
        return


    def get_da_row(self, hit_data, max_n_points):
        hit_lengths = [len(hit_list[1]) for hit_list in hit_data]
        da_row = len(hit_lengths)-1
        for i, length in enumerate(hit_lengths):
            if length < max_n_points:
                da_row = i-1
                break
        if da_row < 0:
            da_row = 0
        return da_row

    def get_title(self, variables):
        title = ""
        for key in variables:
            title += utilities.sub_to_name(key)+": "+str(variables[key])
        return title

    def plot_one_phi(self, row_data, da_key, axis1, axis2, max_n_points, variables):
        global ROOT_OBJECTS
        hit_data = row_data[da_key]
        name = "phi_"+da_key+"_"+axis1+" vs "+axis2
        title = get_title(variables)
        name += " "+title
        canvas = None
        tune_data = []
        area_data = []
        da_row = get_da_row(hit_data, max_n_points)
        print("Plot on phi")
        cholesky_canvas = common.make_root_canvas(name+"_cholesky")
        cholesky_canvas.SetCanvasSize(1000, 1000)
        delta = (da_row+1)*0.4+1.
        hist = common.make_root_histogram("", [-100.], "u", 100, [-100.], "u'", 100,
                                        xmin=-delta, xmax=+delta, ymin=-delta, ymax=+delta)
        hist.Draw()
        for i, hit_list in enumerate(hit_data[:da_row+1]):
            hit_list = [Hit.new_from_dict(hit_dict) for hit_dict in hit_list[1]]
            finder = DPhiTuneFinder()
            finder.u = [hit[axis1] for hit in hit_list[1:]]
            finder.up = [hit[axis2] for hit in hit_list[1:]]
            area_src = numpy.array([[hit[axis1], hit[axis2]] for hit in hit_list])
            try:
                finder.get_tune()
                an_area = get_area(area_src)
            except Exception:
                print("Failed to unpack data for phi da plot")
                continue
            cholesky_canvas, hist, graph = finder.plot_cholesky_space(cholesky_canvas, 1.+i*0.2)
            if i < da_row:
                color = ROOT.TColor(10000+len(ROOT_OBJECTS), 0., 1.-i*1./da_row, 0.)
                ROOT_OBJECTS.append(color)
                graph.SetMarkerColor(ROOT.kGreen)
            ROOT_OBJECTS += [hist, graph]
            tune_data += copy.deepcopy(finder.dphi)
            area_data += [an_area for i in finder.dphi]
            print("    Area", an_area, end=' ')
            print(" tune", numpy.mean(tune_data), "+/-", numpy.std(tune_data))
        canvas = common.make_root_canvas(name)
        canvas.Draw()
        hist, graph = common.make_root_graph(name, area_data, "Amplitude [mm]", tune_data, "Fractional tune")
        ROOT_OBJECTS += [hist, graph]
        hist.SetTitle(da_key)
        hist.Draw()
        graph.SetMarkerStyle(24)
        graph.Draw("p same")
        return canvas, cholesky_canvas

    def plot_acceptance(self, acceptance, x_data, y_data, axis2):
        area_src = numpy.transpose(numpy.array([x_data, y_data]))
        area = self.get_area(area_src)
        scale = (acceptance/area)**0.5
        #covariance *= scale
        x_mean, y_mean = numpy.mean(x_data), numpy.mean(y_data)
        x_values = [(x-x_mean)*scale for x in x_data]
        y_values = [(y-y_mean)*scale for y in y_data]
        print("Acceptance bounding box (relative to mean) x:", min(x_values), max(x_values), "x':", min(y_values), max(y_values))
        values = list(zip(x_values, y_values))
        values = sorted(values, key = lambda x: math.atan2(x[0], x[1]))
        x_values, y_values = list(zip(*values))
        x_values = [x+x_mean for x in x_values]
        y_values = [y+y_mean for y in y_values]
        x_values.append(x_values[0])
        y_values.append(y_values[0])
        hist, graph = common.make_root_graph("", x_values, "", y_values, "", sort = False)
        graph.SetLineColor(ROOT.kOrange)
        return graph

    def plot_one_da(self, row_data, da_key, axis1, axis2, max_n_points, variables, acceptance):
        hit_data = row_data[da_key]
        name = da_key+"_"+axis1+" vs "+axis2
        title = self.get_title(variables)
        name += " "+title
        canvas = None
        da_row = self.get_da_row(hit_data, max_n_points)
        graph_list = []
        acceptance_graph = None
        units = {'x':'mm', 'y':'mm', 'px':'MeV/c', 'py':'MeV/c', "x'":'rad', "y'":'rad'}
        for i, hit_list in enumerate(hit_data):
            hit_list = [Hit.new_from_dict(hit_dict) for hit_dict in hit_list[1]]
            x_data = [hit[axis1] for hit in hit_list]
            y_data = [hit[axis2] for hit in hit_list]
            axis1_label = axis1+" ["+units[axis1]+"]"
            axis2_label = axis2+" ["+units[axis2]+"]"
            hist, graph = common.make_root_graph(name+" "+str(i), x_data, axis1_label, y_data, axis2_label)
            graph.SetMarkerStyle(7)
            if i == da_row:
                canvas = common.make_root_canvas(name)
                hist.SetTitle(title)
                hist.Draw()
            elif i < da_row:
                graph.SetMarkerColor(ROOT.kGreen)
            else:
                graph.SetMarkerColor(ROOT.kGray)
            graph_list.append(graph)
            pz = hit_list[0]['pz']
            mass = common.pdg_pid_to_mass[2212]
            if axis2 == "px" or axis2 == "py":
                geom_acceptance = acceptance*mass #mm mrad normalised -> MeV/c
            elif axis2 == "x'" or axis2 == "y'":
                geom_acceptance = acceptance*mass/pz 
            #if i == 3:
            #    acceptance_graph = self.plot_acceptance(geom_acceptance, x_data, y_data, axis2)

        for graph in graph_list:
            graph.Draw("SAMEP")
        if acceptance_graph != None:
            acceptance_graph.Draw("SAME C")
        canvas.Update()
        return canvas

    def get_seed(self, hit, ref):
        seed = [hit[var]-ref[var] for var in self.var_list]
        return seed

    def make_plot(self, data_list, da_row, name, x_name, y_name):
        graph_list = []
        for x_data, y_data in data_list:
            i = len(graph_list)
            hist, graph = common.make_root_graph(name+str(i), x_data, x_name, y_data, y_name)
            graph.SetMarkerStyle(7)
            if i == da_row:
                canvas = common.make_root_canvas(name)
                hist.SetTitle(name)
                hist.Draw()
            elif i == 0 and i < da_row:
                graph.SetMarkerColor(ROOT.kBlue)
            elif i < da_row:
                graph.SetMarkerColor(ROOT.kGreen)
            else:
                graph.SetMarkerColor(ROOT.kGray)
            graph_list.append(graph)
        canvas.cd()
        canvas.Update()
        for graph in graph_list[da_row+1:]+graph_list[1:da_row]+\
                     graph_list[0:1]+graph_list[da_row:da_row+1]:
            graph.Draw("SAMEP")
        name = name.replace(":", "")
        name = name.replace(" ", "_")
        for fmt in self.formats:
            canvas.Print(self.plot_dir+"/"+name+"."+fmt)

    def plot_decoupled_da(self, row_data, da_key, axis1, axis2, max_n_points, variables, acceptance):
        hit_data = row_data[da_key]
        da_row = self.get_da_row(hit_data, max_n_points)
        tm = row_data["ref_track"][0]
        tm = copy.deepcopy(row_data["tm"])
        for i, row in enumerate(tm):
            tm[i] = row[1:5]
        tm = DecoupledTransferMatrix(tm)
        name = da_key+"_"+axis1+" vs "+axis2
        u_data, v_data, au_data, av_data, a4_data = [], [], [], [], []

        closed_orbit = Hit.new_from_dict(row_data["ref_track"][0])
        for i, hit_list in enumerate(hit_data):
            hit_list = [Hit.new_from_dict(hit_dict) for hit_dict in hit_list[1]]
            seed_list = [self.get_seed(hit, closed_orbit) for hit in hit_list]
            seed_list = [tm.decoupled(seed) for seed in seed_list]
            u = [seed[0] for seed in seed_list]
            pu = [seed[1] for seed in seed_list]
            v = [seed[2] for seed in seed_list]
            pv = [seed[3] for seed in seed_list]
            u_data.append((u, pu))
            v_data.append((v, pv))

            #turns = range(len(seed_list))
            #amp_data_u = [self.get_amplitude(tm, seed, [0, 1]) for seed in seed_list]
            ##amp_data_v = [self.get_amplitude(tm, seed, [2, 3]) for seed in seed_list]
            #amp_data_4d = [self.get_amplitude(tm, seed, [0, 1, 2, 3]) for seed in seed_list]
        self.make_plot(u_data, da_row, da_key+": u vs pu", "u", "pu")
        self.make_plot(v_data, da_row, da_key+": v vs pv", "v", "pv")

    def test(self):
        points = numpy.array([
        [1., -1.], [1., 1.], [-1., 1.], [-1., -1.]
        ])
        area = self.get_area(points)
        print("Square Test - get area should be 4:", area)

        points = numpy.array([
        [2., 0.], [0., 2.], [0., -2.], [-2., 0.], [1., -1.], [1., 1.], [-1., 1.], [-1., -1.]
        ])
        area = self.get_area(points)
        print("Star Test - get area should be 8:", area)

        points += [10., 7.]
        area = self.get_area(points)
        print("Translated Star Test - get area should be 8:", area)

def main():
    base_dir = "output/"
    #base_fname = base_dir+"/baseline/get_da.tmp"
    min_n_turns = 100
    for file_name in glob.glob(sys.argv[1]):
        plotter = DAPlotter(file_name)
        #if file_name != base_fname:
        #    data += load_data(base_fname)
        plot_dir = os.path.split(file_name)[0]
        # acceptance should be normalised mm mrad
        plotter.plot_da(min_n_turns, plot_dir, 2.7*1e-3*math.pi)
        plotter.make_global_plot(min_n_turns, plot_dir) # 500 turn da

if __name__ == "__main__":
    #test()
    main()
    print("Finished - press <CR> to close")
    input()

