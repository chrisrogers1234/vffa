import os
import numpy
import ROOT
import shutil

import xboa.common
import xboa.hit


from opal_tracking import OpalTracking

def clear_dir(dir_name):
    try:
        shutil.rmtree(dir_name)
    except OSError:
        pass
    os.makedirs(dir_name)

def sub_to_name(sub_key):
    sub_name = sub_key[2:-2]
    sub_name = sub_name.replace("_", " ")
    return sub_name


def sub_to_units(sub_key):
    units = {
        "energy":"MeV"
    }
    name = sub_to_name(sub_key)
    if name in units:
        return " ["+units[name]+"]"
    return ""

def reference(config, energy, x=0., px=0.):
    """
    Generate a reference particle
    """
    hit_dict = {}
    hit_dict["pid"] = config.tracking["pdg_pid"]
    hit_dict["mass"] = xboa.common.pdg_pid_to_mass[abs(hit_dict["pid"])]
    hit_dict["charge"] = 1
    hit_dict["x"] = x
    hit_dict["px"] = px
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

def get_substitutions_axis(data):
    subs_ref = data[0]['substitutions']
    axis_candidates = {}
    for item in data:
        subs = item['substitutions']
        for key in subs.keys():
            if subs[key] != subs_ref[key]:
                try:
                    float(subs[key])
                    axis_candidates[key] = []
                except TypeError:
                    continue
    if axis_candidates == {}:
        print "All of", len(data), "items look the same - nothing to plot"
        print "First:"
        print " ", data[0]['substitutions']
        print "Last:"
        print " ", data[-1]['substitutions']
    for item in data:
        for key in axis_candidates:
            axis_candidates[key].append(item['substitutions'][key])
    return axis_candidates

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

