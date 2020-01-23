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
    for line in fin.readlines():
        data_tmp = json.loads(line)
        if data_out == None or data_tmp['score'] < data_out['score']:
            data_out = data_tmp
        n_iterations += 1
    if len(data_out) == 0:
        print(" empty file, abort")
        sys.exit(1)
    print("... loaded", n_iterations, "lines")
    data_out['n_iterations'] = n_iterations
    return [data_out]

def plot_actions(data, plot_dir, foil_probe, tm):
    pass

def plot_score(data, plot_dir, foil_probe, foil_var):
    """
    data - the data
    co_axis - "x" or "y"
    plot_dir - directory to which files are written
    foil_probe - element in the "tracking output" list which contains the foil probe
    """
    print("\nclosed orbit with foil probe", foil_probe)
    metres = 1e-3
    bump_position = []
    bump_score = []
    bump_n_iterations = []
    for bump in data:
        bump_position.append(bump["tracking"][foil_probe][foil_var])
        bump_score.append(bump["score"])
        bump_n_iterations.append(bump["n_iterations"])
    canvas = common.make_root_canvas("score")
    mgraph = ROOT.TMultiGraph()
    root_objects.append(mgraph)
    hist, graph = common.make_root_graph("score",
                                          bump_position*2, var[foil_var],
                                          bump_score+bump_n_iterations, "Score (red)   Iterations (blue)")
    hist, graph = common.make_root_graph("score",
                                          bump_position, "Position [mm]",
                                          bump_score, "Score")
    graph.SetMarkerStyle(20)
    graph.SetMarkerColor(ROOT.kRed)
    mgraph.Add(graph)
    hist, graph = common.make_root_graph("score",
                                          bump_position, var[foil_var],
                                          bump_n_iterations, "n_iterations")
    graph.SetMarkerColor(ROOT.kBlue)
    graph.SetMarkerStyle(24)
    mgraph.Add(graph)
    mgraph.Draw("alp")
    canvas.SetLogy()
    canvas.Update()
    for format in ["png"]:
        canvas.Print(plot_dir+"/score."+format)


def plot_fields(data, plot_dir, foil_probe, foil_var):
    print("\nfield", end=' ')
    bump_fields = {}
    bump_position = []
    for bump_run in data:
        bump_position += [bump["target_bump"][foil_var-1] for bump in bump_run["bumps"]]
        print(bump_position)
        for bump in bump_run["bumps"]:
            for key in bump["bump_fields"]:
                if key not in bump_fields:
                    bump_fields[key] = []
                bump_fields[key].append(bump["bump_fields"][key])
    canvas = common.make_root_canvas("fields")
    all_fields = []
    for key in bump_fields:
        all_fields += bump_fields[key]
    delta = 0.2
    max_pos = (1.+delta)*max(bump_position)-delta*min(bump_position)
    hist, graph = common.make_root_graph("fields",
                                        bump_position*len(bump_fields)+[max_pos],
                                        var[foil_var],
                                        all_fields+[0.], "Dipole field [T]")
    hist.Draw()
    color = [ROOT.kBlue, ROOT.kGreen+2, ROOT.kRed, ROOT.kYellow+2]
    marker = [20, 24, 25, 26]
    leg_list = []
    index = 0
    marker_index = 0
    for key, field in sorted(bump_fields.items()):
        name = key.replace("__", "")
        name = name.replace("_", " ")
        hist, graph = common.make_root_graph(name,
                                            bump_position, var[foil_var],
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


#        score = [bump["score"] for bump in bump_run["bumps"]]
#        n_iterations = [bump["n_iterations"] for bump in bump_run["bumps"]]

def plot_closed_orbit(data, co_axis, plot_dir, station):
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
    pos, mom = [], []
    canvas = common.make_root_canvas("closed orbit "+co_axis+" "+str(station))
    mgraph = ROOT.TMultiGraph()
    marker_styles = [20, 24, 21, 25, 26, 22, 27]
    for i, bump_run in enumerate(data):
        pos, mom = [], []
        for bump in bump_run["bumps"]:
            for item in bump["tracking"]:
                if item[0] == station:
                    pos.append(item[ip])
                    mom.append(item[im])
        hist, graph = common.make_root_graph("closed orbit",
                                             pos, var[ip],
                                             mom, var[im])
        style = marker_styles[i % len(marker_styles)]
        graph.SetMarkerStyle(style)
        mgraph.Add(graph)
    mgraph.Draw("alp")
    root_objects.append(mgraph)
    canvas.Update()
    for format in ["png"]:
        canvas.Print(plot_dir+"/"+co_axis+"_closed_orbit_probe_"+str(station)+"."+format)

root_objects = []
var = {0:"station",
       1:"Radial position [mm]", 2:"Radial momentum [MeV/c]",
       3:"Height [mm]", 4:"Vertical momentum [MeV/c]"}

def main(file_glob):
    foil_probe = 4
    foil_var = 3
    data = []
    for file_name in glob.glob(file_glob):
        print(file_name)
        plot_dir = os.path.split(file_name)[0]+"/plot_bump"
        data += load_file_alt(file_name)

        #for probe in [0, 7]:
        #    for axis in ['x', 'y']:
        #        plot_closed_orbit(data, axis, plot_dir, probe)
        #plot_fields(data, plot_dir, foil_probe, foil_var)
    utilities.clear_dir(plot_dir)
    plot_score(data, plot_dir, foil_probe, foil_var)
    

if __name__ == "__main__":
    main(sys.argv[1])
    input("Done - Press <CR> to end")
