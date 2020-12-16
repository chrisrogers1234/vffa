import glob
import os
import json
import sys

import xboa.common as common
import matplotlib
import matplotlib.pyplot
import ROOT

import utils.utilities as utilities


def load_file_alt(file_name, flip_vertical):
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
    if flip_vertical:
        for hit in data_out["tracking"]:
            hit[3] *= -1
            hit[4] *= -1
    data_out['n_iterations'] = n_iterations
    return [data_out]



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

def plot_all_fields_2d(data, plot_dir, foil_probe, var_1, var_2):
    field_name_list = sorted(data[0]["bump_fields"])
    canvas = ROOT.TCanvas("fields 2d", "fields 2d", 2000, 1000)
    canvas.Divide(4, 2)
    index = 1
    for name in field_name_list:
        this_canvas = canvas.cd(index)
        z_lambda = lambda bump: bump["bump_fields"][name]
        increment = plot_2d(data, foil_probe, 1, 3, z_lambda, name, this_canvas)
        if increment != False:
            index += 1
    root_objects.append(canvas)
    canvas.Print(plot_dir+"/fields_2d.png")

def plot_phased_var_2d(data, plot_dir, foil_probe, var_1, var_2, var_3, number_of_phases):
    canvas = ROOT.TCanvas(var_3, var_3, 2000, 1000)
    canvas.Divide(number_of_phases, 1)
    for index in range(number_of_phases):
        this_canvas = canvas.cd(index+1)
        z_lambda = lambda bump: bump[var_3][index]
        this_canvas.SetLogz()
        this_canvas.Update()
        name = var_3+" for phase "+str(index+1)
        increment = plot_2d(data, foil_probe, 1, 3, z_lambda, name, this_canvas)
    root_objects.append(canvas)
    canvas.Print(plot_dir+"/"+var_3+"_2d.png")

def plot_2d(data, foil_probe, foil_var1, foil_var2, z_lambda, plot_title, canvas):
    name = plot_title.replace("__", "")
    name = name.replace("_", " ")
    x_data, y_data, z_data = [], [], []
    for bump in data:
        for hit in bump["tracking"]:
            if hit[0] == foil_probe:
                x_data.append(hit[foil_var1])
                y_data.append(hit[foil_var2])
                break
        if len(x_data) != len(z_data): # meaning we found foil probe
            z_data.append(z_lambda(bump))
    if abs(max(z_data)-min(z_data)) < 1e-5:
        print("Not plotting", name)
        print("    ", "z_data", z_data)
        return False

    print("Plotting ", name)
    print("   ", "x_data", x_data)
    print("   ", "y_data", y_data)
    print("   ", "z_data", z_data)

    n_points = len(x_data)
    if canvas == None:
        canvas = common.make_root_canvas(name)
        canvas.Draw()
    canvas.cd()
    x_min_max = common.min_max(x_data)
    y_min_max = common.min_max(y_data)
    z_min_max = common.min_max(z_data)
    hist = ROOT.TH3D("", "",
                     10, x_min_max[0], x_min_max[1],
                     10, y_min_max[0], y_min_max[1],
                     10, z_min_max[0], z_min_max[1],
            )
    hist.GetXaxis().SetTitleOffset(3)
    hist.SetStats(False)
    hist.GetZaxis().SetLabelSize(0.01)
    hist.Draw()
    root_objects.append(hist)
    graph2d = ROOT.TGraph2D(n_points)
    root_objects.append(graph2d)
    ROOT.gStyle.SetPalette(1)
    graph2d.SetMarkerStyle(20)
    graph2d.SetTitle("")
    graph2d.Draw('same pcolz')
    for i in range(n_points):
        graph2d.SetPoint(i, x_data[i], y_data[i], z_data[i])
    #palette = ROOT.TPaletteAxis(hist.GetListOfFunctions().FindObject("palette"))
    #palette.SetX2NDC(0.99)
    canvas.SetRightMargin(0.13)
    text = ROOT.TPaveText(0.1, 0.6, 0.65, 0.7)
    text.SetBorderSize(1)
    text.SetFillColor(10)
    text.AddText(name)
    text.Draw()
    root_objects.append(text)
    name = name.replace(" ", "_")
    canvas.SetTheta(90)
    canvas.SetPhi(-90)
    canvas.Update()
    return canvas, graph2d

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
    if len(add_points) > 0:
        pos = [vec[ip] for vec in add_points]
        mom = [vec[im] for vec in add_points]
        pos_all += pos
        mom_all += mom
        hist, graph = common.make_root_graph("closed orbit",
                                             pos, "",
                                             mom, "", sort=False)
        graph.SetLineColor(ROOT.kGray)
        graph.SetMarkerColor(ROOT.kGray)
        style = marker_styles[0]
        graph.SetMarkerStyle(style)
        mgraph.Add(graph)
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
            style = marker_styles[i+1 % len(marker_styles)]
            graph.SetMarkerStyle(style)
            mgraph.Add(graph)
            pos_all += pos
            mom_all += mom
        except IndexError:
            print("Failed to plot")
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

def plot_fields_2(data, plot_dir, x_values):
    bump_fields = {}
    for bump in data:
        for key in bump["bump_fields"]:
            if key not in bump_fields:
                bump_fields[key] = []
            bump_fields[key].append(bump["bump_fields"][key])
    figure = matplotlib.pyplot.figure()
    axes = figure.add_subplot(1, 1, 1)
    if type(x_values) == type(""):
        x_name = x_values.replace("__", "").replace("_", " ").replace("field", "")
        axes.set_xlabel(x_name+" [T]")
        x_values = bump_fields[x_values]
    else:
        axes.set_xlabel("Setting number")
    axes.set_ylabel("Bump field [T]")
    print(bump_fields.keys())
    for key in sorted(bump_fields.keys()):
        name = key.replace("__", "").replace("_", " ").replace("field", "")
        if "h bump" in name:
            style = "dotted"
        else:
            style = "dashed"
        axes.plot(x_values, bump_fields[key], label=name, linestyle=style, marker="o")
    axis_xrange = list(axes.get_xlim())
    axis_xrange[1] += (axis_xrange[1]-axis_xrange[0])*0.5
    axes.set_xlim(axis_xrange)
    axes.legend()
    figure.savefig(plot_dir+"/bump_fields.png")

def load_trajectory(file_name, station, fix, fixed_point, momentum):
    if file_name == None:
        return []
    fin = open(file_name)
    for line in fin.readlines():
        pass
    trajectory = json.loads(line)["beam_trajectory"]
    # offset everything so that trajectory[fixed_point] = fix
    delta = [fix[i] - trajectory[fixed_point][i] for i in range(4)]
    for j, point in enumerate(trajectory):
        point[1] *= momentum
        point[3] *= momentum
        for i in range(4):
            point[i] += delta[i]
        trajectory[j] = [station]+point
    trajectory = trajectory[0:fixed_point]
    return trajectory


root_objects = []
var = {0:"station",
       1:"Radial position [mm]", 2:"Radial momentum [MeV/c]",
       3:"Height [mm]", 4:"Vertical momentum [MeV/c]"}

def main(file_list):
    flip_vertical = False
    number_of_phases = 1
    score_tolerance = 1000.
    foil_probe = 0
    foil_var = 3
    foil_axis = "Height at foil [mm]"
    data = []
    plot_dir = os.path.split(file_list[0])[0]+"/plot_bump"
    utilities.clear_dir(plot_dir)
    index = 0
    score, n_iterations = [], []
    for file_glob in file_list:
        for file_name in sorted(glob.glob(file_glob)):
            index += 1
            if index % number_of_phases != 0:
                score.append(load_file_alt(file_name, flip_vertical)[0]["score"])
                n_iterations.append(load_file_alt(file_name, flip_vertical)[0]["n_iterations"])
                continue
            data += load_file_alt(file_name, flip_vertical)
            data[-1]["score"] = score+[data[-1]["score"]]
            data[-1]["n_iterations"] = n_iterations+[data[-1]["n_iterations"]]
            score, n_iterations = [], []
    traj_fname = "output/triplet_baseline/anticorrelated_painting/toy_model/dp_p_0.0013__col_dens_2e-05__n_tot_50__n_i_25__inj_emit_0.1__k_1.0__tgt_emit_8.0/run_summary.json"
    trajectory = [] #load_trajectory(traj_fname, 0, [3744.07267-10, 3.843, -89.83+8, -1.158], 32, 75.0)
    for probe in [0, 7]:
        for axis in ['x', 'y', 'x-y']:
            if probe == 0:
                add_points = trajectory + \
                             [[probe]+ [3744.07267, 3.843, -89.83, -1.158],] #3744.08063, 3.843, -9.81, 1.151],)
            else:
                add_points = ([probe]+[3739.427804804206, -0.0002874136417290174, -88.77890374233482, -0.0008126002995965109],)
            if flip_vertical:
                for point in add_points:
                    point[2] *= -1
                    point[3] *= -1
            plot_closed_orbit(data, axis, plot_dir, probe, add_points)
    plot_fields(data, plot_dir, foil_probe, foil_var, foil_axis)
    plot_fields_2(data, plot_dir, [0, 1, 2, 3, 4])
    plot_all_fields_2d(data, plot_dir, foil_probe, 3, 1)
    for bump in data:
        print(bump["score"])

    plot_phased_var_2d(data, plot_dir, foil_probe, 3, 1, "score", number_of_phases)
    plot_phased_var_2d(data, plot_dir, foil_probe, 3, 1, "n_iterations", number_of_phases)
    

if __name__ == "__main__":
    main(sys.argv[1:])
    input("Done - Press <CR> to end")

