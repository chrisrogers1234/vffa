"""
Plot a single closed orbit (once tracking has finished)
"""

import os
import sys
import copy
import math
import numpy
from utils import utilities
import plot_dump_fields

MASS = 938.2720813

try:
    import ROOT
except ImportError:
    print "You need to install PyROOT to run this example."

class RootObjects:
    histograms = []
    canvases = []
    graphs = []
    other = []

class Colors:
    ref_colours = numpy.array([ROOT.kRed, ROOT.kOrange+2, ROOT.kYellow+2, 
                          ROOT.kGreen+1, ROOT.kBlue, ROOT.kMagenta+1])
    colours = copy.deepcopy(ref_colours)

    @classmethod
    def next(cls):
        a_colour = cls.colours[0]
        cls.colours = numpy.roll(cls.colours, 1)
        return a_colour
    
    @classmethod
    def reset(cls):
        cls.colours = copy.deepcopy(cls.ref_colours)

def r_phi_track_file(data):
    data = copy.deepcopy(data)
    data["r"] = range(len(data["x"]))
    data["phi"] = range(len(data["x"]))
    data["pr"] = range(len(data["x"]))
    data["pphi"] = range(len(data["x"]))
    for i in range(len(data["r"])):
        data["r"][i] = (data["x"][i]**2+data["y"][i]**2.)**0.5
        phi = math.atan2(data["y"][i], data["x"][i])
        data["phi"][i] = math.degrees(phi)
        if data["phi"][i] < 0.:
            data["phi"][i] += 360.
        p = (data["px"][i]**2+data["py"][i]**2)**0.5
        data["pr"][i] = p*math.sin(phi)
        data["pphi"][i] = p*math.cos(phi)
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
        if "#" in words:
            line = fin.readline()[:-1]
            continue
        if len(words) != len(heading):
            print "Line\n  "+line+"\nmismatched to heading\n  "+str(heading)+"\nin parse_file "+file_name
        else:
            words = [types[i](x) for i, x in enumerate(words)]
            for i, item in enumerate(heading):
                data[item].append(words[i])
        line = fin.readline()[:-1]
    print "Got data from file "+file_name
    return data

def parse_track_file(filename):
    file_name = filename
    heading = ["id", "x", "px", "y", "py", "z", "pz"]#, "bx", "by", "bz", "ex", "ey", "ez"]
    types = [str]+[float]*(len(heading)-1)
    data = parse_file(file_name, heading, types)
    data["px"] = [px*MASS for px in data["px"]]
    data["py"] = [py*MASS for py in data["py"]]
        
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

def plot_x_y_projection(step_list, canvas = None):
    axes = None
    if canvas == None:
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
    graph = ROOT.TGraph(len(step_list))
    for i in range(len(step_list["x"])):
        graph.SetPoint(i, step_list["x"][i], step_list["y"][i])
    graph.Draw("l")
    canvas.Update()
    RootObjects.canvases.append(canvas)
    RootObjects.graphs.append(graph)
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
        RootObjects.histograms.append(axes)
    else:
        canvas.cd()
    graph = ROOT.TGraph(len(step_list))
    points = zip(step_list["phi"], step_list["r"])
    points = sorted(points)
    for i in range(len(step_list["r"])):
        graph.SetPoint(i, points[i][0], points[i][1])
    graph.SetMarkerColor(Colors.next())
    graph.Draw("p")
    canvas.Update()
    RootObjects.canvases.append(canvas)
    RootObjects.graphs.append(graph)
    return canvas, axes, graph

def plot_x_z_projection(step_list, canvas = None):
    if canvas == None:
        canvas = ROOT.TCanvas("x_z_projection", "x_z_projection")
        axes = ROOT.TH2D("x_z_projection_axes", ";#phi [degree];z [m]",
                        1000, 0, 360./5,
                        1000, 0.0, 0.2)
        axes.SetStats(False)
        canvas.Draw()
        axes.Draw()
        RootObjects.histograms.append(axes)
        RootObjects.canvases.append(canvas)
    graph_list = []
    points = zip(step_list["phi"], step_list["z"])
    old_phi = max(step_list["phi"])+1.0
    index = 0
    for i in range(len(step_list["z"])):
        phi, z = points[i]
        if phi < old_phi:
            index = 0
            graph_list.append(ROOT.TGraph())
        graph_list[-1].SetPoint(index, points[i][0], points[i][1])
        index += 1
        old_phi = phi
    for graph in graph_list:
        graph.Draw("l")
    canvas.Update()
    RootObjects.graphs.append(graph)
    return canvas, graph_list

def step_statistics(step_list):
    delta_r_list = []
    for i in range(len(step_list["x"])-1):
        delta_x = step_list["x"][i+1]-step_list["x"][i]
        delta_y = step_list["y"][i+1]-step_list["y"][i]
        delta_z = step_list["z"][i+1]-step_list["z"][i]
        delta_r_list.append((delta_x**2+delta_y**2+delta_z**2)**0.5)
    print len(step_list), "steps with mean size:", numpy.mean(delta_r_list),
    print "and RMS:", numpy.std(delta_r_list)

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

    canvas = None
    for step_list in step_list_of_lists:
        canvas, graph = plot_x_z_projection(step_list, canvas)
    for format in ["png"]:
        canvas.Print(output_dir+"closed_orbit_elevation."+format)
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
    field_plot = plot_dump_fields.PlotDumpFields(opal_run_dir+"FieldMapXY.dat")
    field_plot.load_dump_fields()
    #inner_radius, axis_radius, outer_radius, ncells = 13.5, 14.350, 15.5, 30
    inner_radius, axis_radius, outer_radius, ncells = 3.5, 3.995, 4.5, 20

    canvas = field_plot.plot_dump_fields("x", "y", "bz")
    canvas, axes, graph = plot_x_y_projection(step_list, canvas)
    plot_beam_pipe(inner_radius, outer_radius, ncells, canvas)
    plot_axis(axis_radius, ncells, canvas)
    #plot_elements_xy(opal_run_dir+"log", canvas)
    for format in ["png"]:
        canvas.Print(output_dir+"closed_orbit_plan_bz."+format)

    return
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

    canvas = field_plot.plot_dump_fields("x", "y", "br")
    canvas, axes, graph = plot_x_y_projection(step_list, canvas)
    plot_beam_pipe(inner_radius, outer_radius, ncells, canvas)
    #plot_elements_xy(opal_run_dir+"log", canvas)
    for format in ["png"]:
        canvas.Print(output_dir+"closed_orbit_cartesian_br."+format)

    canvas = field_plot.plot_dump_fields("x", "y", "bphi")
    canvas, axes, graph = plot_x_y_projection(step_list, canvas)
    plot_beam_pipe(inner_radius, outer_radius, ncells, canvas)
    #plot_elements_xy(opal_run_dir+"log", canvas)
    for format in ["png"]:
        canvas.Print(output_dir+"closed_orbit_cartesian_bphi."+format)

    canvas, graph = plot_x_z_projection(step_list)
    for format in ["png"]:
        canvas.Print(output_dir+"closed_orbit_cartesian_vertical."+format)

def plot_test_field(output_dir, opal_run_dir):
    field_plot = plot_dump_fields.PlotDumpFields(opal_run_dir+"FieldMapTest.dat")
    field_plot.load_dump_fields()
    canvas_1d = None
    for field, color in ("bx", 4), ("bz", 1), ("by", 2):
        canvas_1d, graph = field_plot.plot_1d({}, "x", field,
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
                print "Fitting", eqn
                graph.Fit(fit)
                fit.Draw("SAME")
                RootObjects.other.append(fit)
    canvas_1d.Print(output_dir+"test_field.png")

def print_track(tgt_phi, step_list_of_lists):
    for step_list in step_list_of_lists:
        for i, phi in enumerate(step_list['phi']):
            if phi > tgt_phi:
                break
        print "step list item", i
        for key in sorted(step_list):
            print "    ", key, step_list[key][i]
        for j in 0, i:
            print "p_tot at", j, ":",
            p_tot = (step_list['px'][j]**2+step_list['py'][j]**2+step_list['pz'][j]**2)**0.5
            print format(p_tot, "8.4g")
        print

def main(output_dir, run_dir, run_file_list):
    output_dir += "/"
    opal_run_dir = output_dir+run_dir+"/"
    print "OPAL RUN DIR", opal_run_dir
    step_list_of_lists = []
    for run_file in run_file_list:
        step_list_of_lists.append(parse_track_file(opal_run_dir+run_file))
    plot_cartesian(output_dir, opal_run_dir, step_list_of_lists[0])
    #plot_cylindrical(output_dir, opal_run_dir, step_list_of_lists)
    #plot_test_field(output_dir, opal_run_dir)
    #print_track(0.1*360./15, step_list_of_lists)
    return

    step_statistics(step_list)

if __name__ == "__main__":
    utilities.setup_gstyle()
    tgt_dir = os.path.split(sys.argv[1])[0]
    output_dir = os.path.split(tgt_dir)[0]
    run_dir = os.path.split(tgt_dir)[1]
    run_file_list = [os.path.split(arg)[1] for arg in sys.argv[1:]]
    main(output_dir, run_dir, run_file_list)

