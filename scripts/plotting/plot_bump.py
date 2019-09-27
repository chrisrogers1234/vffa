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
    for bump_run in data:
        bump_position += [bump["tracking"][foil_probe][foil_var] for bump in bump_run["bumps"]]
        bump_score += [bump["score"] for bump in bump_run["bumps"]]
        bump_n_iterations += [bump["n_iterations"] for bump in bump_run["bumps"]]
    canvas = common.make_root_canvas("score")
    hist, graph = common.make_root_graph("score",
                                          bump_position*2, var[foil_var],
                                          bump_score+bump_n_iterations, "Score (red)   Iterations (blue)")
    hist.Draw()
    hist, graph = common.make_root_graph("score",
                                          bump_position, "Radial position [mm]",
                                          bump_score, "Score")
    graph.SetMarkerStyle(20)
    graph.SetMarkerColor(ROOT.kRed)
    graph.Draw("lpsame")
    hist, graph = common.make_root_graph("score",
                                          bump_position, var[foil_var],
                                          bump_n_iterations, "n_iterations")
    graph.SetMarkerColor(ROOT.kBlue)
    graph.SetMarkerStyle(24)
    graph.Draw("lpsame")
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

def plot_closed_orbit(data, co_axis, plot_dir, foil_probe):
    """
    data - the data
    co_axis - "x" or "y"
    plot_dir - directory to which files are written
    foil_probe - element in the "tracking output" list which contains the foil probe
    """
    print("\nclosed orbit with foil probe", foil_probe)
    if co_axis == "x":
        ip = 1 # position index
        im = 2 # momentum index
    elif co_axis == "y":
        ip = 3
        im = 4
    pos, mom = [], []
    for bump_run in data:
        pos += [bump["tracking"][foil_probe][ip] for bump in bump_run["bumps"]]
        mom += [bump["tracking"][foil_probe][im] for bump in bump_run["bumps"]]
        canvas = common.make_root_canvas("closed orbit")
        hist, graph = common.make_root_graph("closed orbit",
                                             pos, var[ip],
                                             mom, var[im])
        hist.Draw()
    marker_styles = [20, 24, 21, 25, 26, 22, 27]
    for i, bump_run in enumerate(data):
        pos = [bump["tracking"][foil_probe][ip] for bump in bump_run["bumps"]]
        mom = [bump["tracking"][foil_probe][im] for bump in bump_run["bumps"]]
        hist, graph = common.make_root_graph("closed orbit",
                                             pos, var[ip],
                                             mom, var[im])
        style = marker_styles[i % len(marker_styles)]
        graph.SetMarkerStyle(style)
        graph.Draw("lpsame")
    canvas.Update()
    for format in ["png"]:
        canvas.Print(plot_dir+"/"+co_axis+"_closed_orbit_probe_"+str(foil_probe)+"."+format)

var = {0:"station",
       1:"Radial position [mm]", 2:"Radial momentum [MeV/c]",
       3:"Height [mm]", 4:"Vertical momentum [MeV/c]"}

def main(file_glob):
    foil_probe = 4
    foil_var = 3
    for file_name in glob.glob(file_glob):
        plot_dir = os.path.split(file_name)[0]
        data = load_file(file_name)
        for probe in range(foil_probe-2, foil_probe+4):
            for axis in ['x', 'y']:
                plot_closed_orbit(data, axis, plot_dir, probe)
        plot_fields(data, plot_dir, foil_probe, foil_var)
        plot_score(data, plot_dir, foil_probe, foil_var)
    

if __name__ == "__main__":
    main(sys.argv[1])
    print("Done")
