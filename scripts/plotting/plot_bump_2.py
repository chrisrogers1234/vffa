import glob
import os
import json
import sys

import xboa.common as common
import ROOT

import utils.utilities as utilities



def load_file(file_name):
    fin = open(file_name)
    print("Loading", file_name, end=' ')
    sys.stdout.flush()
    data_out = json.loads(fin.read())
    if len(data_out) == 0:
        print(" empty file, abort")
        sys.exit(1)
    print(" done")
    return data_out

def load_file_alt(file_name):
    fin = open(file_name)
    print("Loading", file_name, end=' ')
    sys.stdout.flush()
    data_out = None
    n_iterations = 0
    score_list = []
    for line in fin.readlines():
        data_tmp = json.loads(line)
        score_list.append(data_tmp['score'])
        if data_out == None or data_tmp['score'] < data_out['score']:
            data_out = data_tmp
        n_iterations += 1
    if len(data_out) == 0:
        print(" empty file, abort")
        sys.exit(1)
    print("... loaded", n_iterations, "lines with",
          "worst score", format(max(score_list), "8.4g"), 
          "best score", format(data_out['score'], "8.4g"))
    data_out['n_iterations'] = n_iterations
    return [data_out]

def plot_actions(data, plot_dir, foil_probe, tm):
    pass

def get_bump_position(data, foil_probe, foil_var):
    bump_position = []
    for bump in data:
        bump_position.append(0.)
        for hit in bump["tracking"]:
            if hit[0] == foil_probe:
                bump_position[-1] = hit[foil_var]
                break
    return bump_position

def plot_score(data, plot_dir, foil_probe, foil_var, axis_name):
    """
    data - the data
    co_axis - "x" or "y"
    plot_dir - directory to which files are written
    foil_probe - element in the "tracking output" list which contains the foil probe
    """
    print("\nclosed orbit with foil probe", foil_probe)
    metres = 1e-3
    bump_score = []
    bump_n_iterations = []
    bump_position = get_bump_position(data, foil_probe, foil_var)
    for bump in data:
        bump_score.append(bump["score"])
        bump_n_iterations.append(bump["n_iterations"])
    canvas = common.make_root_canvas("score")
    mgraph = ROOT.TMultiGraph()
    root_objects.append(mgraph)
    hist, graph = common.make_root_graph("score",
                                          bump_position*2, axis_name,
                                          bump_score+bump_n_iterations, "Score (red)   Iterations (blue)")
    hist, graph = common.make_root_graph("score",
                                          bump_position, "Position [mm]",
                                          bump_score, "Score")
    graph.SetMarkerStyle(20)
    graph.SetMarkerColor(ROOT.kRed)
    mgraph.Add(graph)
    mgraph.GetXaxis().SetTitle(axis_name)
    hist, graph = common.make_root_graph("Number of iterations",
                                          bump_position, axis_name,
                                          bump_n_iterations, "n_iterations")
    graph.SetMarkerColor(ROOT.kBlue)
    graph.SetMarkerStyle(24)
    mgraph.Add(graph)
    mgraph.Draw("alp")
    #canvas.BuildLegend()
    canvas.SetLogy()
    canvas.Update()
    for format in ["png"]:
        canvas.Print(plot_dir+"/score."+format)


def plot_fields(data, plot_dir, foil_probe, foil_var, axis_name):
    bump_fields = {}
    bump_position = get_bump_position(data, foil_probe, foil_var)
    for bump in data:
        for key in bump["bump_fields"]:
            if key not in bump_fields:
                bump_fields[key] = []
            bump_fields[key].append(bump["bump_fields"][key])
    print("bump_fields = {")
    for key in bump_fields:
        print("    \""+key+"\" :", bump_fields[key], ",")
    print("}")
    canvas = common.make_root_canvas("fields")
    all_fields = []
    for key in bump_fields:
        all_fields += bump_fields[key]
    delta = 0.2
    max_pos = (1.+delta)*max(bump_position)-delta*min(bump_position)
    hist, graph = common.make_root_graph("fields",
                                        bump_position*len(bump_fields)+[max_pos],
                                        axis_name,
                                        all_fields+[0.], "Dipole field [T]")
    hist.Draw()
    color = [ROOT.kBlue, ROOT.kGreen+2, ROOT.kRed, ROOT.kYellow+2, ROOT.kCyan+1]
    marker = [20, 24, 21, 25]
    leg_list = []
    index = 0
    marker_index = 0
    for key, field in sorted(bump_fields.items()):
        name = key.replace("__", "")
        name = name.replace("_", " ")
        hist, graph = common.make_root_graph(name,
                                            bump_position, axis_name,
                                            field, "Dipole field [T]")
        graph.SetMarkerStyle(marker[marker_index])
        graph.SetMarkerColor(color[index])
        graph.SetLineColor(color[index])
        graph.SetLineStyle(color[index])
        graph.Draw("lpsame")
        leg_list.append(graph)
        index += 1
        if index == len(color):
            index = 0
            marker_index += 1
    leg = common.make_root_legend(canvas, leg_list)
    leg.SetBorderSize(1)
    leg.SetX1NDC(0.75)
    leg.SetX2NDC(0.9)
    leg.SetY1NDC(0.5)
    leg.SetY2NDC(0.9)
    leg.Draw()
    canvas.Update()
    for format in ["png"]:
        canvas.Print(plot_dir+"/bump_fields."+format)

def get_limits(data_list, tolerance=1e-9):
    min_x = min(data_list)
    max_x = max(data_list)
    if max_x - min_x < tolerance:
        max_x += 0.5
        min_x -= 0.5
    else:
        delta = max_x-min_x
        max_x += delta/10.
        min_x -= delta/10.
    return min_x, max_x

def plot_closed_orbit(data, co_axis, plot_dir, station, add_points=()):
    """
    data - the data
    co_axis - "x" or "y"
    plot_dir - directory to which files are written
    station - station in the "tracking output" list at which to plot
    """
    print("\nclosed orbit with foil probe", station)
    if co_axis == "x":
        ip = 1 # position index
        im = 2 # momentum index
    elif co_axis == "y":
        ip = 3
        im = 4
    elif co_axis == "x-y":
        ip = 1
        im = 3
    pos_all, mom_all = [], []
    canvas = common.make_root_canvas("closed orbit "+co_axis+" "+str(station))
    mgraph = ROOT.TMultiGraph()
    marker_styles = [20, 20, 21, 21, 22, 24, 24, 25, 25, 26]
    for i, bump in enumerate(data):
        pos, mom = [], []
        station_list = [item[0] for item in bump["tracking"]]
        try:
            index = station_list.index(station)
            pos.append(bump["tracking"][index][ip])
            mom.append(bump["tracking"][index][im])
        except ValueError:
            print("Failed to find station", station, "in list", station_list)
        try:
            hist, graph = common.make_root_graph("closed orbit",
                                                    pos, "",
                                                    mom, "")
            style = marker_styles[i % len(marker_styles)]
            graph.SetMarkerStyle(style)
            mgraph.Add(graph)
            pos_all += pos
            mom_all += mom
        except IndexError:
            print("Failed to plot")
    if len(add_points) > 0:
        pos = [vec[ip] for vec in add_points]
        mom = [vec[im] for vec in add_points]
        pos_all += pos
        mom_all += mom
        hist, graph = common.make_root_graph("closed orbit",
                                             pos, "",
                                             mom, "")
        graph.SetMarkerColor(ROOT.kBlue)
        style = marker_styles[i % len(marker_styles)]
        graph.SetMarkerStyle(style)
        mgraph.Add(graph)
    mgraph.GetXaxis().SetTitle(var[ip])
    mgraph.GetYaxis().SetTitle(var[im])
    pos_limits = get_limits(pos_all, 0.1)
    mom_limits = get_limits(mom_all, 0.1)
    mgraph.GetXaxis().SetLimits(pos_limits[0], pos_limits[1]);
    mgraph.SetMinimum(mom_limits[0]);
    mgraph.SetMaximum(mom_limits[1]);
    mgraph.Draw("alp")
    root_objects.append(mgraph)
    canvas.Update()
    for format in ["png"]:
        canvas.Print(plot_dir+"/"+co_axis+"_closed_orbit_probe_"+str(station)+"."+format)

root_objects = []
var = {0:"station",
       1:"Radial position [mm]", 2:"Radial momentum [MeV/c]",
       3:"Height [mm]", 4:"Vertical momentum [MeV/c]"}

def main(file_list):
    foil_probe = 0
    foil_var = 3
    foil_axis = "Height at foil [mm]"
    data = []
    plot_dir = os.path.split(file_list[0])[0]+"/plot_bump"
    utilities.clear_dir(plot_dir)
    index = 0
    for file_glob in file_list:
        for file_name in glob.glob(file_glob):
            index += 1
            if index %2 == 1:
                continue
            print(file_name)
            data += load_file_alt(file_name)
    for probe in [0, 5]:
        for axis in ['x', 'y', 'x-y']:
            if probe == 0:
                add_points = ([0, 3869.70412, -19.79, 137.8, 3.752],)
            else:
                add_points = ([6, 3946.122, -24.619, 136.662, -1.236],)
            plot_closed_orbit(data, axis, plot_dir, probe, add_points)
    plot_fields(data, plot_dir, foil_probe, foil_var, foil_axis)
    plot_score(data, plot_dir, foil_probe, foil_var, foil_axis)
    

if __name__ == "__main__":
    main(sys.argv[1:])
    input("Done - Press <CR> to end")
