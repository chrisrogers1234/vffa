import glob
import os
import json
import sys

import xboa.common as common
import ROOT

import utilities



def load_file(file_name):
    fin = open(file_name)
    print "Loading", file_name,
    sys.stdout.flush()
    data_out = json.loads(fin.read())
    if len(data_out) == 0:
        print " empty file, abort"
        sys.exit(1)
    print " done"
    return data_out

def plot_score(data, plot_dir, foil_probe):
    """
    data - the data
    co_axis - "x" or "y"
    plot_dir - directory to which files are written
    foil_probe - element in the "tracking output" list which contains the foil probe
    """
    print "\nclosed orbit with foil probe", foil_probe
    metres = 1e-3
    bump_position = []
    bump_score = []
    bump_n_iterations = []
    for bump_run in data:
        bump_position += [bump["tracking"][foil_probe][1] for bump in bump_run["bumps"]]
        bump_score += [bump["score"] for bump in bump_run["bumps"]]
        bump_n_iterations += [bump["n_iterations"] for bump in bump_run["bumps"]]
    canvas = common.make_root_canvas("score")
    hist, graph = common.make_root_graph("score",
                                          bump_position*2, "Radial position [mm]",
                                          bump_score+bump_n_iterations, "Score (red)   Iterations (blue)")
    hist.Draw()
    hist, graph = common.make_root_graph("score",
                                          bump_position, "Radial position [mm]",
                                          bump_score, "Score")
    graph.SetMarkerStyle(20)
    graph.SetMarkerColor(ROOT.kRed)
    graph.Draw("lpsame")
    hist, graph = common.make_root_graph("score",
                                          bump_position, "Radial position [mm]",
                                          bump_n_iterations, "n_iterations")
    graph.SetMarkerColor(ROOT.kBlue)
    graph.SetMarkerStyle(24)
    graph.Draw("lpsame")
    canvas.SetLogy()
    canvas.Update()
    for format in "eps", "png", "root":
        canvas.Print(plot_dir+"/score."+format)


def plot_fields(data, plot_dir, foil_probe):
    print "\nfield",
    bump_fields = [[], [], [], []]
    bump_position = []
    for bump_run in data:
        bump_position += [bump["tracking"][foil_probe][1] for bump in bump_run["bumps"]]
        for bump in bump_run["bumps"]:
            for i in range(4):
                bump_fields[i].append(bump["bump_fields"][i])
    canvas = common.make_root_canvas("fields")
    all_fields = bump_fields[0]+bump_fields[1]+bump_fields[2]+bump_fields[3]
    delta = 0.2
    max_pos = (1.+delta)*max(bump_position)-delta*min(bump_position)
    hist, graph = common.make_root_graph("fields",
                                        bump_position*4+[max_pos], "Radial position [mm]",
                                        all_fields+[0.], "Dipole field [T]")
    hist.Draw()
    color = [ROOT.kBlue, ROOT.kGreen+2, ROOT.kRed, ROOT.kYellow+2]
    leg_list = []
    for i, field in enumerate(bump_fields):
        hist, graph = common.make_root_graph("Bump "+str(i+1),
                                            bump_position, "Radial position [mm]",
                                            field, "Dipole field [T]")
        graph.SetMarkerStyle(20)
        graph.SetMarkerColor(color[i])
        graph.SetLineColor(color[i])
        graph.Draw("lpsame")
        leg_list.append(graph)
    leg = common.make_root_legend(canvas, leg_list)
    leg.SetBorderSize(1)
    leg.SetX1NDC(0.75)
    leg.SetX2NDC(0.9)
    leg.SetY1NDC(0.5)
    leg.SetY2NDC(0.9)
    leg.Draw()
    canvas.Update()
    for format in "eps", "png", "root":
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
    print "\nclosed orbit with foil probe", foil_probe
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
                                             pos, "Radial position [mm]",
                                             mom, "Radial momentum [MeV/c]")
        hist.Draw()
    marker_styles = [20, 24, 21, 25, 26, 22, 27]
    for i, bump_run in enumerate(data):
        pos = [bump["tracking"][foil_probe][ip] for bump in bump_run["bumps"]]
        mom = [bump["tracking"][foil_probe][im] for bump in bump_run["bumps"]]
        hist, graph = common.make_root_graph("closed orbit",
                                             pos, "Radial position [mm]",
                                             mom, "Radial momentum [MeV/c]")
        style = marker_styles[i % len(marker_styles)]
        graph.SetMarkerStyle(style)
        graph.Draw("lpsame")
    canvas.Update()
    for format in "eps", "png", "root":
        canvas.Print(plot_dir+"/"+co_axis+"_closed_orbit_probe_"+str(foil_probe)+"."+format)

def main():
    foil_probe = 5
    file_glob = "output/bump_design_no_bump/find_bump_parameters.tmp"
    for file_name in glob.glob(file_glob):
        plot_dir = os.path.split(file_name)[0]
        data = load_file(file_name)
        for probe in range(foil_probe-2, foil_probe+4):
            plot_closed_orbit(data, 'x', plot_dir, probe)
        plot_fields(data, plot_dir, foil_probe)
        plot_score(data, plot_dir, foil_probe)
    

if __name__ == "__main__":
    main()
    print "Done"
