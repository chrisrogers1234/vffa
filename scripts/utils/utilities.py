import copy
import os
import shutil
import io
import sys
import importlib

import numpy
import ROOT

import xboa.common
import xboa.hit


from opal_tracking import OpalTracking
from opal_tracking import StoreDataInMemory
from opal_tracking import PyOpalTracking

import matplotlib
import matplotlib.pyplot

def sub_to_name(sub_key):
    sub_name = sub_key[2:-2]
    sub_name = sub_name.replace("_", " ")
    sub_name = sub_name[0].upper()+sub_name[1:]
    return sub_name

def sub_to_units(sub_key):
    units = {
        "energy":"MeV"
    }
    name = sub_to_name(sub_key).lower()
    if name in units:
        return " ["+units[name]+"]"
    return ""

def clear_dir(a_dir):
    try:
        shutil.rmtree(a_dir)
    except OSError:
        pass
    os.makedirs(a_dir)

def preprocess_subs(subs):
    subs_tmp = copy.deepcopy(subs)
    for key, value in subs_tmp.items():
        if type(value) == type([]):
            list_str = "{"
            for list_item in value[:-1]:
                list_str += str(list_item)+", "
            list_str += str(value[-1])+"}"
            subs[key] = list_str
        elif value is True:
            subs[key] = "TRUE"
        elif value is False:
            subs[key] = "FALSE"
    return subs


def do_lattice(config, subs, overrides):
    subs = copy.deepcopy(subs)
    subs.update(overrides)
    subs = preprocess_subs(subs)
    lattice_in = config.tracking["lattice_file"]
    lattice_out = config.tracking["lattice_file_out"]
    xboa.common.substitute(lattice_in, lattice_out, subs)
    return subs

def reference(config, energy, x=0., px=0., y=0., py=0.):
    """
    Generate a reference particle
    """
    hit_dict = {}
    hit_dict["pid"] = config.tracking["pdg_pid"]
    hit_dict["mass"] = xboa.common.pdg_pid_to_mass[abs(hit_dict["pid"])]
    hit_dict["charge"] = 1
    hit_dict["x"] = x
    hit_dict["px"] = px
    hit_dict["y"] = y
    hit_dict["py"] = py
    hit_dict["kinetic_energy"] = energy
    hit = xboa.hit.Hit.new_from_dict(hit_dict, "pz")
    return hit

def setup_tracking(config, probes, ref_energy):
    ref_hit = reference(config, ref_energy)
    opal_exe = os.path.expandvars(config.tracking["opal_path"])
    lattice = config.tracking["lattice_file_out"]
    log = config.tracking["tracking_log"]
    beam = config.tracking["beam_file_out"]
    tracking = OpalTracking(lattice, beam, ref_hit, probes, opal_exe, log)
    tracking.verbose = config.tracking["verbose"]
    tracking.set_file_format(config.tracking["file_format"])
    tracking.flags = config.tracking["flags"]
    tracking.pass_through_analysis = StoreDataInMemory(config)
    return tracking

PY_OPAL_TRACKING = None
def setup_py_tracking(config, run_dir, phi_list):
    global PY_OPAL_TRACKING
    if PY_OPAL_TRACKING is not None:
        print("Warning - tried to setup PyOpal again - abort")
        return PY_OPAL_TRACKING
    lattice = config.tracking["lattice_file_out"]
    tracking = PyOpalTracking(config, run_dir)
    tracking.step_list = phi_list
    PY_OPAL_TRACKING = tracking
    return tracking

def tune_lines(canvas, min_order=0, max_order=8):
    canvas.cd()
    for x_power in range(0, max_order):
        for y_power in range(0, max_order):
            if y_power + x_power > max_order or y_power + x_power == 0:
                continue
            x_points = [0., 1.]
            y_points = [0., 1.]
            if x_power > y_power:
                x_points[0] = y_points[0]*y_power/x_power
                x_points[1] = y_points[1]*y_power/x_power
            else:
                y_points[0] = x_points[0]*x_power/y_power
                y_points[1] = x_points[1]*x_power/y_power
            hist, graph = xboa.common.make_root_graph("", x_points, "", y_points, "")
            graph.Draw("SAMEL")
            x_points = [1.-x for x in x_points]
            #y_points = [1.-y for y in y_points]
            hist, graph = xboa.common.make_root_graph("", x_points, "", y_points, "")
            graph.Draw("SAMEL")
    canvas.Update()

def get_substitutions_axis(data, subs_key):
    subs_ref = data[0][subs_key]
    axis_candidates = {}
    for item in data:
        subs = item[subs_key]
        for key in list(subs.keys()):
            try:
                comp = subs[key] == subs_ref[key]
            except KeyError:
                print("Warning - missing ", key, "treating as same")
                comp = True
            if not comp:
                try:
                    float(subs[key])
                    axis_candidates[key] = []
                except (TypeError, ValueError):
                    continue
    if axis_candidates == {}:
        print("All of", len(data), "items look the same - nothing to plot")
        print("First:")
        print(" ", data[0][subs_key])
        print("Last:")
        print(" ", data[-1][subs_key])
    for item in data:
        for key in axis_candidates:
            #print key, axis_candidates.keys()
            #print item.keys()
            #print item[subs_key].keys()
            axis_candidates[key].append(item[subs_key][key])
    return axis_candidates

def get_groups(data, group_axis, subs_key):
    axis_candidates = get_substitutions_axis(data, subs_key)
    if group_axis != None and group_axis not in axis_candidates:
        raise RuntimeError("Did not recognise group axis "+str(group_axis))
    if group_axis == None:
        group_list = [item for item in data]
        return {"":{"item_list":[i for i in range(len(data))]}}
    else:
        # group_list is lists the possible groups
        # e.g. list of all possible values of "__bump_field_1__"
        group_list = [item[subs_key][group_axis] for item in data]
        group_list = list(set(group_list)) # unique list
        # tmp_group_dict is mapping from group value to the items having that value
        tmp_group_dict = dict([(group, []) for group in group_list])
        for key in tmp_group_dict:
            tmp_group_dict[key] = [i for i, item in enumerate(data) \
                                        if item[subs_key][group_axis] == key]
    group_dict = {}
    for key in tmp_group_dict:
        new_key = sub_to_name(group_axis)+" "+format(key, "3.3g")
        group_dict[new_key] = {'item_list':tmp_group_dict[key]}
        print(new_key, ":", group_dict[new_key])
    return group_dict

def matplot_marker_size(x_points):
    marker_size = 1
    if len(x_points):
        marker_size = 100/len(x_points)**0.5
    return marker_size

def setup_gstyle():
    stops = [0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000]
    red   = [0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764]
    green = [0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832]
    blue  = [0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539]
    s = numpy.array(stops)
    r = numpy.array(red)
    g = numpy.array(green)
    b = numpy.array(blue)

    ncontours = 255
    npoints = len(s)
    ROOT.TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
    ROOT.gStyle.SetNumberContours(ncontours)

    # axes and labels
    ROOT.gStyle.SetPadBottomMargin(0.15)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    ROOT.gStyle.SetPadRightMargin(0.15)
    for axis in "X", "Y":
        ROOT.gStyle.SetNdivisions(505, axis)
        ROOT.gStyle.SetLabelSize(0.05, axis)
        ROOT.gStyle.SetTitleSize(0.06, axis)
        ROOT.gStyle.SetTitleOffset(1.10, axis)

def setup_da_figure(include_projections):
    fig = matplotlib.pyplot.figure(figsize=(20, 10))
    y, dy = 0.10, 0.35
    if include_projections:
        y, dy = 0.15, 0.30
    axes = [
        fig.add_subplot(2, 3, 1,  position=[0.06, 0.55, 0.26, 0.35]),
        fig.add_subplot(2, 6, 7,  position=[0.06, y, 0.10, dy]),
        fig.add_subplot(2, 6, 8,  position=[0.22, y, 0.10, dy]),
        fig.add_subplot(2, 3, 2,  position=[0.38, 0.55, 0.26, 0.35]),
        fig.add_subplot(2, 6, 9,  position=[0.38, 0.10, 0.10, 0.35]),
        fig.add_subplot(2, 6, 10, position=[0.54, 0.10, 0.10, 0.35]),
        fig.add_subplot(2, 3, 3,  position=[0.70, 0.55, 0.26, 0.35]),
        fig.add_subplot(2, 6, 11, position=[0.70, 0.10, 0.10, 0.35]),
        fig.add_subplot(2, 6, 12, position=[0.86, 0.10, 0.10, 0.35]),
    ]

    if include_projections:
        axes += [
            fig.add_subplot(2, 7, 13,  position=[0.06, 0.10, 0.10, 0.05]),
            fig.add_subplot(2, 7, 14,  position=[0.22, 0.10, 0.10, 0.05]),
        ]
    return fig, axes

def plot_id(figure, lattice, study, version, x=0.9, y=0.95):
    text = lattice.replace("_", " ")
    text += "\n"+study.replace("_", " ")
    if version != "release":
        text += " "+str(version)
    else:
        text += " "+str(version)
    figure.text(x, y, text, fontsize=10)

def directory_name(lattice, study, version):
    dir_name = lattice.replace(" ", "_")+"_"+study.replace(" ", "_")
    if version == "release":
        dir_name += "r"
    else:
        dir_name += "_"+str(version)
    dir_name += "/"
    return dir_name

def get_config():
    if len(sys.argv) < 2:
        print("Usage: python /path/to/run_one.py /path/to/config.py")
        sys.exit(1)
    config_file_name = sys.argv[1]
    config_file = sys.argv[1].replace(".py", "")
    config_file = config_file.split("scripts/")[1]
    config_file = config_file.replace("/", ".")
    config_args = tuple(sys.argv[2:])
    print("Using configuration module", config_file, "with arguments", config_args)
    config_mod = importlib.import_module(config_file)
    config = config_mod.Config(*config_args)
    return config_file_name, config
