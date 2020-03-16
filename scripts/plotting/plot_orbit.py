"""
Plot a single closed orbit (once tracking has finished)
"""

import argparse
import os
import sys
import copy
import math
import glob
import numpy
import h5py

from utils import utilities
import plotting.plot_dump_fields as plot_dump_fields
import PyOpal.parser
import PyOpal.field


MASS = 938.2720813
TARGET = "ID1"

try:
    import ROOT
except ImportError:
    print("You need to install PyROOT to run this example.")

class RootObjects:
    histograms = []
    canvases = []
    graphs = []
    other = []

class Colors:
    ref_colours = numpy.array([ROOT.kBlue, ROOT.kGreen+1, ROOT.kRed, 
                               ROOT.kOrange+2, ROOT.kYellow+2, ROOT.kMagenta+1])
    colours = copy.deepcopy(ref_colours)

    @classmethod
    def next(cls):
        a_colour = cls.colours[0]
        cls.colours = numpy.roll(cls.colours, -1)
        return int(a_colour)

    @classmethod
    def current(cls):
        return int(cls.colours[0])

    @classmethod
    def reset(cls):
        cls.colours = copy.deepcopy(cls.ref_colours)

class CoordinateTransform:
    def __init__(self, cell_angle, cell_length):
        self.cell_angle = cell_angle
        self.cell_length = cell_length
        self.r0 = cell_length/2./math.sin(cell_angle/2.)
        self.phi0 = math.pi/2.

    def phi(self, x, y):
        return math.atan2(-x, y)

    def r(self, x, y):
        return (x**2+y**2)**0.5

    def cell_number(self, x, y):
        return math.floor(self.phi(x, y)/self.cell_angle)

    def cell_phi(self, x, y):
        return self.cell_number(x, y)*self.cell_angle+self.cell_angle/2.

    def cell_start(self, x, y):
        phi_start = self.cell_number(x, y)*self.cell_angle
        y0 = self.r0*math.cos(phi_start)
        x0 = -self.r0*math.sin(phi_start)
        return x0, y0

    def coordinate_transform(self, x, y, bx, by, bz):
        x0, y0 = self.cell_start(x, y)
        cell_phi = self.cell_phi(x, y)
        cell_number = self.cell_number(x, y)
        cell_x = cell_number*self.cell_length

        x_cell = -math.cos(cell_phi)*(x-x0)-math.sin(cell_phi)*(y-y0)+cell_x
        y_cell = -math.sin(cell_phi)*(x-x0)+math.cos(cell_phi)*(y-y0)
        bx_cell = +math.cos(cell_phi)*bx+math.sin(cell_phi)*by
        by_cell = +math.sin(cell_phi)*bx-math.cos(cell_phi)*by
        bz_cell = -bz
        #print(x, y, "**", math.degrees(cell_phi), "**", x0, y0, "**", x_cell, y_cell)
        return x_cell, y_cell, bx_cell, by_cell, bz_cell

    @classmethod
    def test_coordinate_transform(cls):
        transform = CoordinateTransform(math.radians(12), 3.0)
        x, y, bx, by, bz = transform.coordinate_transform(-1e-9, transform.r0, 1., 0., 1.)
        print(x, y, bx, by, bz)
        assert(abs(x) < 1e-6)
        assert(abs(y) < 1e-6)
        assert(abs(-bz-1.0) < 1e-6)
        assert(abs(bx-math.cos(transform.cell_angle/2.)) < 1e-6)
        assert(abs(by-math.sin(transform.cell_angle/2.)) < 1e-6)
        x, y, bx, by, bz = transform.coordinate_transform(-1e-12, transform.r0*1e-6, 1., 0., 1.)

def r_phi_track_file(data):
    data = copy.deepcopy(data)
    data["r"] = list(range(len(data["x"])))
    data["phi"] = list(range(len(data["x"])))
    data["pr"] = list(range(len(data["x"])))
    data["pphi"] = list(range(len(data["x"])))
    for i in range(len(data["r"])):
        data["r"][i] = (data["x"][i]**2+data["y"][i]**2.)**0.5
        phi = math.atan2(data["y"][i], data["x"][i])
        data["phi"][i] = math.degrees(phi)
        if data["phi"][i] < 0.:
            data["phi"][i] += 360.
        px = data["px"][i]
        py = data["py"][i]
        #p = (data["px"][i]**2+data["py"][i]**2)**0.5
        data["pr"][i]   = px*math.cos(phi)+py*math.sin(phi)
        data["pphi"][i] = -px*math.sin(phi)+py*math.cos(phi)
    return data

def parse_file(file_name, heading, types):
    if len(heading) != len(types):
        raise KeyError("Heading mismatched to types in parse_file "+file_name)
    fin = open(file_name)
    data = {}
    for item in heading:
        data[item] = []
    line = fin.readline()[:-1]
    while line != "":
        words = line.split()
        if  "#" in words or (TARGET != None and words[0] != TARGET):
            line = fin.readline()[:-1]
            continue
        if len(words) != len(heading):
            print("Line\n  "+line+"\nmismatched to heading\n  "+str(heading)+"\nin parse_file "+file_name)
        else:
            words = [types[i](x) for i, x in enumerate(words)]
            for i, item in enumerate(heading):
                data[item].append(words[i])
        line = fin.readline()[:-1]
    print("Got data from file "+file_name)
    return data

def parse_track_file(filename):
    file_name = filename
    heading = ["id", "x", "px", "y", "py", "z", "pz"]#, "bx", "by", "bz", "ex", "ey", "ez"]
    types = [str]+[float]*(len(heading)-1)
    data = parse_file(file_name, heading, types)
    data["px"] = [px*MASS for px in data["px"]]
    data["py"] = [py*MASS for py in data["py"]]
    data["pz"] = [pz*MASS for pz in data["pz"]]
    data = r_phi_track_file(data)
    return data

def load_track_orbit(file_name):
    fin = open(file_name)
    step_list = []
    for i, line in enumerate(fin.readlines()):
        if i < 2:
          continue
        words = line.split()
        step = {}
        step["particle_index"] = int(words[0][2:])
        step["x_pos"] = float(words[1])
        step["beta_x_gamma"] = float(words[2])
        step["y_pos"] = float(words[3])
        step["beta_y_gamma"] = float(words[4])
        step["z_pos"] = float(words[5])
        step["beta_z_gamma"] = float(words[6])
        step_list.append(step)
    return step_list

def plot_x_y_projection(step_list_of_lists, canvas = None):
    axes = None
    if canvas == None:
        Colors.reset()
        canvas = ROOT.TCanvas("x_y_projection", "x_y_projection")
        canvas.Draw()
        axes = ROOT.TH2D("x_y_projection_axes", ";x [m];y [m]",
                         1000, -25., 25.,
                         1000, -25., 25.)
        axes.SetStats(False)
        axes.Draw()
        RootObjects.histograms.append(axes)
    else:
        canvas.cd()
    for step_list in step_list_of_lists:
        graph = ROOT.TGraph(len(step_list))
        for i in range(len(step_list["x"])):
            graph.SetPoint(i, step_list["x"][i], step_list["y"][i])
        graph.SetMarkerColor(Colors.next())
        graph.SetLineColor(Colors.current())
        graph.Draw("l")
        RootObjects.graphs.append(graph)
    canvas.Update()
    RootObjects.canvases.append(canvas)
    return canvas, axes, graph

def plot_r_phi_projection(step_list, canvas = None):
    axes = None
    if canvas == None:
        canvas = ROOT.TCanvas("x_y_projection", "x_y_projection")
        canvas.Draw()
        axes = ROOT.TH2D("x_y_projection_axes", ";#phi [degree];r [m]",
                         1000, 90., 150.,
                         1000, 1., 4.)
        axes.SetStats(False)
        axes.Draw()
        Colors.reset()
        RootObjects.histograms.append(axes)
    else:
        canvas.cd()
    graph = ROOT.TGraph(len(step_list))
    points = list(zip(step_list["phi"], step_list["r"]))
    points = sorted(points)
    for i in range(len(step_list["r"])):
        graph.SetPoint(i, points[i][0], points[i][1])
    graph.SetMarkerColor(Colors.next())
    graph.Draw("p")
    canvas.Update()
    RootObjects.canvases.append(canvas)
    RootObjects.graphs.append(graph)
    return canvas, axes, graph

def plot_x_z_projection(step_list_of_lists, phi0, phi1, z0, z1, canvas = None):
    if canvas == None:
        canvas = ROOT.TCanvas("x_z_projection", "x_z_projection")
        axes = ROOT.TH2D("x_z_projection_axes", ";#phi [degree];z [m]",
                        1000, phi0, phi1,
                        1000, z0, z1)
        axes.SetStats(False)
        canvas.Draw()
        axes.Draw()
        RootObjects.histograms.append(axes)
        RootObjects.canvases.append(canvas)
    graph_list = []
    for step_list in step_list_of_lists:
        graph_list.append(ROOT.TGraph())
        graph_list[-1].SetLineColor(Colors.next())
        points = list(zip(step_list["phi"], step_list["z"]))
        old_phi = max(step_list["phi"])+1.0
        index = 0
        for i in range(len(step_list["z"])):
            phi, z = points[i]
            if phi < old_phi and abs(phi-old_phi) > 90.0:
                index = 0
                graph_list.append(ROOT.TGraph())
                graph_list[-1].SetLineColor(Colors.current())
            graph_list[-1].SetPoint(index, points[i][0], points[i][1])
            index += 1
            old_phi = phi
    for graph in graph_list:
        graph.Draw("l SAME")
        RootObjects.graphs.append(graph)
    canvas.Update()
    return canvas, graph_list

def step_statistics(step_list):
    delta_r_list = []
    for i in range(len(step_list["x"])-1):
        delta_x = step_list["x"][i+1]-step_list["x"][i]
        delta_y = step_list["y"][i+1]-step_list["y"][i]
        delta_z = step_list["z"][i+1]-step_list["z"][i]
        delta_r_list.append((delta_x**2+delta_y**2+delta_z**2)**0.5)
    print(len(step_list), "steps with mean size:", numpy.mean(delta_r_list), end=' ')
    print("and RMS:", numpy.std(delta_r_list))

def plot_phi_pipe(n_periods, phi0, phi1, y0, y1, color, canvas):
    line_style=2
    canvas.cd()
    print("PLOT PHI PIPE", n_periods, y0, y1, phi0, phi1, color)
    for i in range(n_periods):
        phi = phi0+(phi1-phi0)/float(n_periods-1)*i
        graph = ROOT.TGraph(2)
        graph.SetPoint(0, phi, y0)
        graph.SetPoint(1, phi, y1)
        graph.SetLineStyle(line_style)
        graph.SetLineColor(color)
        graph.Draw("l SAME")
        RootObjects.graphs.append(graph)
    graph = ROOT.TGraph(2)
    graph.SetPoint(0, phi0, y0)
    graph.SetPoint(1, phi1, y0)
    graph.SetLineStyle(line_style)
    graph.SetLineColor(color)
    graph.Draw("l SAME")
    RootObjects.graphs.append(graph)
    graph = ROOT.TGraph(2)
    graph.SetPoint(0, phi0, y1)
    graph.SetPoint(1, phi1, y1)
    graph.SetLineStyle(line_style)
    graph.SetLineColor(color)
    graph.Draw("l SAME")
    RootObjects.graphs.append(graph)

def plot_beam_pipe(inner_radius, outer_radius, n_periods, canvas=None):
    n_steps = 361 # number of azimuthal steps

    if canvas == None:
        canvas = ROOT.TCanvas("beam_pipe", "beam_pipe")
        canvas.Draw()
        axes = ROOT.TH2D("beam_pipe_axes", ";x [mm];y [mm]",
                         1000, -25., 25.,
                         1000, -25., 25.)
        axes.Draw()
        RootObjects.histograms.append(axes)
        RootObjects.canvases.append(canvas)
    canvas.cd()
    graph_inner = ROOT.TGraph(n_steps)
    graph_outer = ROOT.TGraph(n_steps)
    for i in range(n_steps):
        graph_inner.SetPoint(i,
                             inner_radius*math.sin(i/float(n_steps-1)*2.*math.pi),
                             inner_radius*math.cos(i/float(n_steps-1)*2.*math.pi))
        graph_outer.SetPoint(i,
                             outer_radius*math.sin(i/float(n_steps-1)*2.*math.pi),
                             outer_radius*math.cos(i/float(n_steps-1)*2.*math.pi))
    for i in range(n_periods):
        graph = ROOT.TGraph(2)
        graph.SetPoint(0,
                       inner_radius*math.sin(i/float(n_periods)*2.*math.pi),
                       inner_radius*math.cos(i/float(n_periods)*2.*math.pi))
        graph.SetPoint(1,
                       outer_radius*math.sin(i/float(n_periods)*2.*math.pi),
                       outer_radius*math.cos(i/float(n_periods)*2.*math.pi))
        graph.Draw("l")
        RootObjects.graphs.append(graph)
    graph_inner.Draw("l")
    graph_outer.Draw("l")
    canvas.Update()
    RootObjects.graphs.append(graph_inner)
    RootObjects.graphs.append(graph_outer)

def plot_axis(axis_radius, n_periods, canvas=None):
    if canvas == None:
        canvas = ROOT.TCanvas("beam_pipe", "beam_pipe")
        canvas.Draw()
        axes = ROOT.TH2D("beam_pipe_axes", ";x [mm];y [mm]",
                         1000, -25., 25.,
                         1000, -25., 25.)
        axes.Draw()
        RootObjects.histograms.append(axes)
        RootObjects.canvases.append(canvas)
    canvas.cd()
    graph_axis = ROOT.TGraph(n_periods+1)
    for i in range(n_periods+1):
        if (n_periods/4)*4 == n_periods:
            index = i
        else:
            index = i+0.5
        phi = index*math.pi*2./n_periods
        graph_axis.SetPoint(i,
                            axis_radius*math.cos(phi),
                            axis_radius*math.sin(phi))

    RootObjects.graphs.append(graph_axis)
    graph_axis.SetLineColor(ROOT.kGray)
    graph_axis.Draw("l")
    canvas.Update()

def load_probes(probe_files, cuts):
    probe_data = []
    for a_file in sorted(glob.glob(probe_files)):
        if ".h5" in a_file:
            probe_data += load_h5_probe(a_file, cuts)
        else:
            probe_data += load_ascii_probe(a_file, cuts)
    return probe_data

def load_ascii_probe(a_file, cuts):
    data = []
    fin = open(a_file)
    fin.readline()
    for line in fin.readlines():
        words = line.split()
        numbers = []
        for a_word in words:
            try:
                a_word = float(a_word)
                numbers.append(a_word)
            except ValueError:
                pass
        if len(numbers) == 0:
            continue
        item = {
            "x":numbers[0],
            "y":numbers[1],
            "z":numbers[2],
            "phi":math.degrees(math.atan2(numbers[1], numbers[0])),
            "r":(numbers[0]**2+numbers[1]**2)**0.5
        }
        will_cut = [cut_lambda(item) for cut_lambda in cuts]
        if sum(will_cut):
            continue
        data.append(item)
    return data

def load_h5_probe(file_name, cuts):
    data = []
    print("Loading h5 file", file_name, end=" ")
    sys.stdout.flush()
    h5_file = h5py.File(file_name, 'r')
    for key in h5_file.keys():
        if key[:5] != "Step#":
            print("skipping key", key, end=" ")
            continue
        n_steps = len(h5_file[key]["x"])
        h5_step = h5_file[key]
        for i in range(n_steps):
            item = {
                "id":h5_step["id"][i],
                "x":h5_step["x"][i],
                "y":h5_step["y"][i],
                "z":h5_step["z"][i],
                "phi":math.degrees(math.atan2(h5_step["y"][i],
                                              h5_step["x"][i])),
                "r":(h5_step["x"][i]**2+h5_step["y"][i]**2)**0.5
            }
            will_cut = [cut_lambda(item) for cut_lambda in cuts]
            if sum(will_cut):
                continue
            data.append(item)
    print("...", len(data), "points")
    return data

def plot_probes(canvas, probe_data, axis_1, axis_2):
    canvas.cd()
    graph = ROOT.TGraph(0)
    point = 0
    for i, item in enumerate(probe_data):
        graph.SetPoint(i, item[axis_1], item[axis_2])
    graph.SetMarkerStyle(24)
    graph.Draw("P SAME")
    graph.SetName("probes_"+axis_1+"_"+axis_2)
    RootObjects.graphs.append(graph)

def plot_b_field(step_list):
    canvas = ROOT.TCanvas("bfield", "bfield")
    axes = ROOT.TH2D("bfield_axes", ";phi [rad];b [T]",
                     1000, -2.*math.pi, 2.*math.pi,
                     1000, -25., 25.)
    axes.SetStats(False)
    graph = ROOT.TGraph(len(step_list))
    canvas.Draw()
    axes.Draw()
    for i, step in enumerate(step_list):
        graph.SetPoint(i, step["x_pos"], step["y_pos"])
    graph.Draw("l")
    canvas.Update()
    RootObjects.histograms.append(axes)
    RootObjects.canvases.append(canvas)
    RootObjects.graphs.append(graph)
    return canvas, axes, graph

def scrape_vector(line):
    vector = line.split("(")[1]
    vector = vector.split(")")[0]
    vector = [float(number)/1000. for number in vector.split(",")]
    return vector

def get_elements(log_file):
    log_file = open(log_file)
    start_positions = []
    for line in log_file.readlines():
        if "Start position (" not in line:
            continue
        start_positions.append(scrape_vector(line))
    return start_positions
  
def plot_elements_xz(log_file, canvas):
    start_positions = get_elements(log_file)
    canvas.cd()
    graph = ROOT.TGraph(len(start_positions))
    for i, pos in enumerate(start_positions):
        phi = math.atan2(pos[1], pos[0])
        graph.SetPoint(i, phi, pos[2])
    RootObjects.graphs.append(graph)
    graph.SetMarkerStyle(7)
    graph.Draw("PSAME")
    canvas.Update()
    return graph

def plot_elements_xy(log_file, canvas):
    start_positions = get_elements(log_file)
    canvas.cd()
    graph = ROOT.TGraph(len(start_positions))
    for i, pos in enumerate(start_positions):
        graph.SetPoint(i, pos[0], pos[1])
    RootObjects.graphs.append(graph)
    graph.SetMarkerStyle(7)
    graph.Draw("PSAME")
    canvas.Update()
    return graph

def plot_cylindrical(output_dir, opal_run_dir, step_list_of_lists):
    field_plot = plot_dump_fields.PlotDumpFields(opal_run_dir+"FieldMapRPHI.dat", True)
    field_plot.load_dump_fields()
    canvas_bz_offset = field_plot.plot_dump_fields("phi", "r", "bz")
    for step_list in step_list_of_lists:
        plot_r_phi_projection(step_list, canvas_bz_offset)
    for format in ["png"]:
        canvas_bz_offset.Print(output_dir+"closed_orbit_cylindrical_bz."+format)
    Colors.reset()

    return
    canvas_br_offset = field_plot.plot_dump_fields("phi", "r", "br")
    for step_list in step_list_of_lists:
        plot_r_phi_projection(step_list, canvas_br_offset)
    for format in ["png"]:
        canvas_br_offset.Print(output_dir+"closed_orbit_cylindrical_br."+format)
    Colors.reset()

    canvas_bphi_offset = field_plot.plot_dump_fields("phi", "r", "bphi")
    for step_list in step_list_of_lists:
        plot_r_phi_projection(step_list, canvas_bphi_offset)
    for format in ["png"]:
        canvas_bphi_offset.Print(output_dir+"closed_orbit_cylindrical_bphi."+format)
    Colors.reset()

    canvas_bphi_offset = field_plot.plot_dump_fields("phi", "r", "bx")
    for step_list in step_list_of_lists:
        plot_r_phi_projection(step_list, canvas_bphi_offset)
    for format in ["png"]:
        canvas_bphi_offset.Print(output_dir+"closed_orbit_cylindrical_bx."+format)
    Colors.reset()

    canvas_bphi_offset = field_plot.plot_dump_fields("phi", "r", "by")
    for step_list in step_list_of_lists:
        plot_r_phi_projection(step_list, canvas_bphi_offset)
    for format in ["png"]:
        canvas_bphi_offset.Print(output_dir+"closed_orbit_cylindrical_by."+format)
    Colors.reset()

    try:
        canvas_1d = None
        for field, color in ("bz", 1), ("bphi", 2), ("br", 4):
            canvas_1d, graph = field_plot.plot_1d({"r":4.}, "phi", field,
                                                  canvas_1d, (0., 36.), color)
        for format in ["png"]:
            canvas_1d.Print(output_dir+"bfield_1d_a."+format)
        canvas_1d = None
        for field, color in ("bz", 1), ("bphi", 2), ("br", 4):
            canvas_1d, graph = field_plot.plot_1d({"r":4.}, "phi", field,
                                                  canvas_1d, (108., 144.), color)
        for format in ["png"]:
            canvas_1d.Print(output_dir+"bfield_1d_b."+format)
        Colors.reset()
    except Exception:
        sys.excepthook(*sys.exc_info())
    


def plot_zoom(output_dir, opal_run_dir, step_list_of_lists):
    field_plot = plot_dump_fields.PlotDumpFields(opal_run_dir+"FieldMapXY-zoom.dat")
    field_plot.load_dump_fields()

    canvas = field_plot.plot_dump_fields("x", "y", "bz")
    for step_list in step_list_of_lists:
        canvas, axes, graph = plot_x_y_projection(step_list, canvas)
    plot_beam_pipe(2.7, 3.7, 2, canvas)
    #plot_elements_xy(opal_run_dir+"log", canvas)
    for format in ["png"]:
        canvas.Print(output_dir+"closed_orbit_plan-zoom."+format)

def plot_cartesian(output_dir, opal_run_dir, step_list):
    field_plot = plot_dump_fields.PlotDumpFields(opal_run_dir+"FieldMapXY.dat", "magnetic_cartesian") #"em_cartesian")
    field_plot.load_dump_fields()

    probe_cuts = [lambda item: item["id"] == 0] # cut if any of probe_cuts is true
    probe_data = load_probes(opal_run_dir+"*PROBE*.h5", probe_cuts)
    #inner_radius, axis_radius, outer_radius, ncells = 12.5, 14.350, 16.0, 30
    #inner_radius_a, axis_radius_a, outer_radius_a, ncells_a = 11.5, 14.350, 16.0, 15
    #z_min, z_max = 0.05, 0.2
    inner_radius, axis_radius, outer_radius, ncells = 3.5, 3.995, 4.5, 20
    inner_radius_a, axis_radius_a, outer_radius_a, ncells_a = 3.4, 3.995, 4.5, 10
    z_min, z_max = -0.0, 0.3

    Colors.reset()
    canvas = field_plot.plot_dump_fields("x", "y", "bz")
    log_plot = plot_dump_fields.LogFileStripper(opal_run_dir+"log",
                                                ["MAGNETD", "MAGNETF"])
    log_plot.plot_log_file(canvas, 0, 1)
    plot_beam_pipe(inner_radius, outer_radius, ncells, canvas)
    plot_beam_pipe(inner_radius_a, outer_radius_a, ncells_a, canvas)
    plot_axis(axis_radius, ncells, canvas)
    plot_probes(canvas, probe_data, "x", "y")
    canvas, axes, graph = plot_x_y_projection(step_list, canvas)
    for format in ["png"]:
        canvas.Print(output_dir+"closed_orbit_plan_bz."+format)
    return
    Colors.reset()
    canvas = None
    canvas, graph = plot_x_z_projection(step_list, 0, 360, z_min, z_max, canvas)
    plot_probes(canvas, probe_data, "phi", "z")
    plot_phi_pipe(ncells+1, 0, 360, -2.0, 2.0, ROOT.kGray, canvas)
    plot_phi_pipe(ncells_a+1, 0, 360, -2.0, 2.0, 1, canvas)
    plot_phi_pipe(2, 72-36*0.07, 72+36*0.07, 0.05, 0.08, ROOT.kGray, canvas)
    plot_phi_pipe(2, 0, 360, 0.0878, -0.05, ROOT.kGray, canvas)
    for format in ["png"]:
        canvas.Print(output_dir+"closed_orbit_elevation."+format)
    Colors.reset()

    Colors.reset()
    canvas = field_plot.plot_dump_fields("x", "y", "br")
    plot_beam_pipe(inner_radius, outer_radius, ncells, canvas)
    plot_beam_pipe(inner_radius_a, outer_radius_a, ncells_a, canvas)
    plot_axis(axis_radius, ncells, canvas)
    plot_probes(canvas, probe_data, "x", "y")
    canvas, axes, graph = plot_x_y_projection(step_list, canvas)
    for format in ["png"]:
        canvas.Print(output_dir+"closed_orbit_plan_br."+format)

    Colors.reset()
    canvas = field_plot.plot_dump_fields("x", "y", "bphi")
    plot_beam_pipe(inner_radius, outer_radius, ncells, canvas)
    plot_beam_pipe(inner_radius_a, outer_radius_a, ncells_a, canvas)
    plot_axis(axis_radius, ncells, canvas)
    plot_probes(canvas, probe_data, "x", "y")
    canvas, axes, graph = plot_x_y_projection(step_list, canvas)
    for format in ["png"]:
        canvas.Print(output_dir+"closed_orbit_plan_bphi."+format)
    return

    Colors.reset()
    canvas = field_plot.plot_dump_fields("x", "y", "btot")
    plot_beam_pipe(inner_radius, outer_radius, ncells, canvas)
    plot_beam_pipe(inner_radius_a, outer_radius_a, ncells_a, canvas)
    plot_axis(axis_radius, ncells, canvas)
    plot_probes(canvas, probe_data, "x", "y")
    canvas, axes, graph = plot_x_y_projection(step_list, canvas)
    for format in ["png"]:
        canvas.Print(output_dir+"closed_orbit_plan_btot."+format)

    canvas = field_plot.plot_dump_fields("x", "y", "bx")
    canvas, axes, graph = plot_x_y_projection(step_list, canvas)
    plot_beam_pipe(inner_radius, outer_radius, ncells, canvas)
    #plot_elements_xy(opal_run_dir+"log", canvas)
    for format in ["png"]:
        canvas.Print(output_dir+"closed_orbit_cartesian_bx."+format)

    canvas = field_plot.plot_dump_fields("x", "y", "by")
    canvas, axes, graph = plot_x_y_projection(step_list, canvas)
    plot_beam_pipe(inner_radius, outer_radius, ncells, canvas)
    #plot_elements_xy(opal_run_dir+"log", canvas)
    for format in ["png"]:
        canvas.Print(output_dir+"closed_orbit_cartesian_by."+format)
    canvas, graph = plot_x_z_projection(step_list)
    for format in ["png"]:
        canvas.Print(output_dir+"closed_orbit_cartesian_vertical."+format)

def plot_test_field(output_dir, opal_run_dir):
    field_plot = plot_dump_fields.PlotDumpFields(opal_run_dir+"FieldMapTest.dat")
    field_plot.load_dump_fields()
    canvas_1d = None
    for field, color in ("bx", 4), ("bz", 1), ("by", 2):
        canvas_1d, graph = field_plot.plot_1d({}, "y", field,
                                                canvas_1d, (-2, 2), color)
        if field == "bz":
            for delta in -1.5, +1.5:
                eqn = "0.5*[0]*(tanh((x-0.5*[1]-[3])/[2])-tanh((x+0.5*[1]-[3])/[2]))"
                fit = ROOT.TF1("fa1", eqn, delta-0.5, delta+0.5)
                fit.SetParameter(0, 0.25)
                fit.SetParameter(1, 0.125)
                fit.SetParameter(2, 0.1)
                fit.SetParameter(3, delta)
                fit.SetLineColor(ROOT.kGreen)
                fit.SetLineStyle(3)
                print("Fitting", eqn)
                graph.Fit(fit)
                fit.Draw("SAME")
                RootObjects.other.append(fit)
    canvas_1d.Print(output_dir+"test_field.png")

def print_track(tgt_phi, step_list_of_lists):
    for step_list in step_list_of_lists:
        for i, phi in enumerate(step_list['phi']):
            if phi > tgt_phi:
                break
        print("step list item", i)
        for key in sorted(step_list):
            print("    ", key, step_list[key][i])
        for j in 0, i:
            print("p_tot at", j, ":", end=' ')
            p_tot = (step_list['px'][j]**2+step_list['py'][j]**2+step_list['pz'][j]**2)**0.5
            print(format(p_tot, "8.4g"))
        print()

def load_lattice(lattice_file):
    if lattice_file == None:
        return
    here = os.getcwd()
    a_dir, a_file = os.path.split(lattice_file)
    os.chdir(a_dir)
    PyOpal.parser.initialise_from_opal_file(a_file)
    os.chdir(here)

def get_machida_field():
    import plot_machida_field
    canvas = plot_machida_field.main()
    return canvas

def plot_orbit_field(output_dir, step_list_of_lists, canvas):
    output = open(output_dir+"/rogers_tracking.txt", "w")
    transform = CoordinateTransform(math.radians(18), 1.25)
    a_list = step_list_of_lists[0]
    n_steps = len(a_list['x'])
    z0 = a_list["z"][0]
    name_list = ["Bx", "By", "Bz"]#, "y", "(z-"+str(z0)+")*100"]
    g_list = [ROOT.TGraph(n_steps) for name in name_list]
    var_list = ['x', 'y', 'z']
    min_b = 0. # minimum axis delta
    min_b_read = 0.
    max_b_read = 0.
    max_b = 20.
    cell_old, x_cell, y_cell = 0., 0., 0.
    x_init = 0.
    print("               x/m               y/m               z/m               bx/T              by/T              bz/T         ", file=output)
    for i in range(n_steps):
        [x, y, z] = [a_list[var][i] for var in var_list]
        try:
            oob, bx, by, bz, dummy, dummy, dummy = \
                                    PyOpal.field.get_field_value(x, y, z, 0.)
        except ValueError:
            print("Could not plot orbit field - no field object available")
            return
        cell_number = transform.cell_number(x, y)
            
        #if abs(cell_number-cell_old) > 1e-9:
        #    print("Cell transform pre-step from", cell_old, "to", cell_number, "step", i)
        #    print("   ", x_cell, y_cell)
        x_cell, y_cell, bx_cell, by_cell, bz_cell = \
                                transform.coordinate_transform(x, y, bx, by, bz)
        if i == 0:
            x_init = x_cell
        x_cell = x_cell-x_init
        #if abs(cell_number-cell_old) > 1e-9:
        #    print("Cell transform post-step")
        #    print("   ", x_cell, y_cell)
        cell_old = cell_number
        if abs(bx) > max_b:
            bx = max_b*bx/abs(bx)
        if abs(by) > max_b:
            by = max_b*by/abs(by)
        if abs(bz) > max_b:
            bz = max_b*bz/abs(bz)
        bm = (bx**2+by**2)**0.5
        g_list[0].SetPoint(i, x_cell, -bx_cell)
        g_list[1].SetPoint(i, x_cell, -by_cell)
        g_list[2].SetPoint(i, x_cell, -bz_cell)
        min_b_read = min([-bx_cell, -by_cell, -bz_cell, min_b_read])
        max_b_read = max([-bx_cell, -by_cell, -bz_cell, max_b_read])
        #g_list[3].SetPoint(i, x_cell, y_cell)
        #g_list[0].SetPoint(i, x_cell, (z-z0))
        print(format(x_cell, "16.8g"),
              format(y_cell, "16.8g"),
              format(z, "16.8g"),
              format(-bx_cell, "16.8g"),
              format(-by_cell, "16.8g"),
              format(-bz_cell, "16.8g"), file=output)
        if i == 0:
            print(x, y, z, "**", bx, by, bz)
    print

    multigraph = ROOT.TMultiGraph()
    multigraph.SetTitle(";s [m];B [T]")
    Colors.reset()
    for i, graph in enumerate(g_list):
        graph.SetTitle(name_list[i])
        graph.SetLineWidth(2)
        graph.SetLineColor(Colors.next())
    if canvas == None:
        canvas = ROOT.TCanvas("Field", "Field")
        canvas.Draw()
        for i, graph in enumerate(g_list):
            multigraph.Add(graph)
        multigraph.Draw("L AXIS")
        if -min_b < min_b_read:
            multigraph.SetMinimum(-min_b)
        if min_b > max_b_read:
            multigraph.SetMaximum(min_b)
        multigraph.GetXaxis().SetNdivisions(510)
        multigraph.GetYaxis().SetNdivisions(510)
    else:
        canvas.cd()
        for i, graph in enumerate(g_list):
            graph.Draw("L SAME")
    RootObjects.graphs.append(multigraph)
    RootObjects.canvases.append(canvas)
    canvas.BuildLegend()
    canvas.Print(output_dir+"/event_field.png")

def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument('file_name_list', type=str, nargs='+')
    parser.add_argument('--lattice_file', dest='lattice_file')
    args = parser.parse_args()
    tgt_dir = os.path.split(args.file_name_list[0])[0]
    output_dir = os.path.split(tgt_dir)[0]
    run_dir = os.path.split(tgt_dir)[1]
    run_file_list = [os.path.split(arg)[1] for arg in args.file_name_list]
    lattice_file = args.lattice_file
    return output_dir, run_dir, run_file_list, lattice_file

def main(output_dir, run_dir, run_file_list, lattice_file):
    load_lattice(lattice_file)
    canvas = None #get_machida_field()
    output_dir += "/"
    opal_run_dir = output_dir+run_dir+"/"
    print("OPAL RUN DIR", opal_run_dir)
    step_list_of_lists = []
    for run_file in run_file_list:
        step_list_of_lists.append(parse_track_file(opal_run_dir+run_file))
    plot_cartesian(output_dir, opal_run_dir, step_list_of_lists)
    #plot_cylindrical(output_dir, opal_run_dir, step_list_of_lists)
    try:
        plot_orbit_field(output_dir, step_list_of_lists, canvas)
    except ValueError:
        sys.excepthook(*sys.exc_info())
    return

    step_statistics(step_list)

if __name__ == "__main__":
    utilities.setup_gstyle()
    args = parse_args(sys.argv)
    #CoordinateTransform.test_coordinate_transform()
    main(*args)
    input("Press <CR> to finish")
