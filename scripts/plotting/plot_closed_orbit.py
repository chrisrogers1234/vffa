import glob
import os
import json
import sys

import xboa.common as common

import utilities

def load_file(file_name):
    fin = open(file_name)
    data_out = []
    for line in fin.readlines():
        data_out.append(json.loads(line))
    return data_out

def plot_closed_orbit(data, co_axis, plot_dir):
    print "\nclosed orbit",
    axis_candidates = utilities.get_substitutions_axis(data)
    metres = 1e-3
    for key in axis_candidates:
        sys.stdout.flush()
        x_name = utilities.sub_to_name(key)
        x_list = axis_candidates[key]
        y_list = [item["hits"][0][co_axis] for item in data]
        y_list = [y*metres for y in y_list]
        canvas = common.make_root_canvas("closed orbit vs "+x_name)
        x_name += utilities.sub_to_units(key)
        hist, graph = common.make_root_graph("closed orbit", x_list, x_name, y_list, "Radial position [m]")
        hist.Draw()
        graph.SetMarkerStyle(4)
        graph.Draw("psame")
        canvas.Update()
        for format in "eps", "png", "root":
            canvas.Print(plot_dir+"/"+co_axis+"_closed_orbit."+format)

def main():
    for file_name in glob.glob("output/rogers_hack/find_closed_orbit.out"):
        plot_dir = os.path.split(file_name)[0]
        print file_name
        data = load_file(file_name)
        if len(data) == 0:
            continue
        plot_closed_orbit(data, 'x', plot_dir) # or ring tune
    

if __name__ == "__main__":
    main()
    raw_input("Done")
