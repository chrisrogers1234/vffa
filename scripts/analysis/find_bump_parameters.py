"""
Script to find the RF set up; drives find closed orbit
"""

import os
import sys
import copy
import json
import math

import numpy
import ROOT

import xboa.common
from xboa.hit import Hit
sys.path.insert(1, "scripts")
from opal_tracking import OpalTracking
import find_closed_orbits
from utils import utilities

class FindBumpParameters(object):
    def __init__(self, config):
        self.config = config
        self.tmp_dir = os.path.join(self.config.run_control["output_dir"],
                               self.config.find_bump_parameters["run_dir"])
        self.x_err = self.config.find_bump_parameters["position_tolerance"]
        self.p_err = self.config.find_bump_parameters["momentum_tolerance"]
        self.first_tracking = True
        self.subs = {}
        self.iteration = 0
        self.target_bump = [0., 0.]
        self.tracking_result = []
        self.output = []
        self.score = None
        self.target_hit = None
        self.overrides = {}

    def setup_minuit(self):
        self.iteration = 0
        self.minuit = ROOT.TMinuit(4)
        min_b = self.config.find_bump_parameters["magnet_min_field"]
        max_b = self.config.find_bump_parameters["magnet_max_field"]
        fixed = self.config.find_bump_parameters["fix_bumps"]
        seed = self.config.find_bump_parameters["seed_field"]
        errors = self.config.find_bump_parameters["seed_errors"]
        for i in range(4):
            parameter_name = "__bump_field_"+str(i+1)+"__"
            self.minuit.DefineParameter(i, parameter_name, seed[i], abs(errors[i]), min_b, max_b)
            if parameter_name in fixed:
                self.minuit.FixParameter(i)
                print "Fixed bump", parameter_name
        self.minuit.SetFCN(self.function)

    def find_bump_parameters(self):
        for subs in self.config.substitution_list:
            self.output.append({"subs":subs, "bumps":[]})
            for self.target_bump in self.config.find_bump_parameters["bump"]:
                self.setup_minuit()
                self.subs = subs
                self.overrides = self.config.find_bump_parameters["subs_overrides"]
                self.max_iterations = self.config.find_bump_parameters["max_iterations"]
                try:
                    self.minuit.Command("SIMPLEX "+str(self.max_iterations)+" 1.")
                except Exception:
                    sys.excepthook(*sys.exc_info())
                    print "Minuit failed"
                self.overrides = self.config.find_bump_parameters["final_subs_overrides"]
                self.track_one(self.get_fields_from_minuit())
                self.save_state("tmp", True)
        self.save_state("out", False)

    def get_filename_root(self):
        fname = self.config.run_control["output_dir"]+"/"
        fname += self.config.find_bump_parameters["output_file"]
        return fname

    def save_state(self, suffix, append):
        state = {
            "target_bump":self.target_bump,
            "bump_fields":self.get_fields_from_minuit(),
            "tracking":self.tracking_result,
            "n_iterations":self.iteration,
            "target_hit":self.target_hit,
            "score":self.score,
        }
        if append:
            self.output[-1]["bumps"].append(state)
        fname = self.get_filename_root()+"."+suffix
        fout = open(fname, "w")
        print >> fout, json.dumps(self.output, indent=2)
 
    def get_fields_from_minuit(self):
        fields = [None]*4
        for i in range(4):
            var = ROOT.Double()
            err = ROOT.Double()
            self.minuit.GetParameter(i, var, err)
            fields[i] = var
        return fields

    def get_score(self, hit_list):
        ignore_stations = self.config.find_bump_parameters["ignore_stations"]
        target_co = self.config.find_bump_parameters["closed_orbit"]
        bump_probe_station = self.config.find_bump_parameters["bump_probe_station"]
        n_hits = len(hit_list)-len(ignore_stations)
        x_score, px_score = 0., 0.
        for i, hit in enumerate(hit_list):
            print format(hit["station"], "6"), format(hit["x"], "8.5g"), format(hit["px"], "8.5g"), format(hit["t"], "12.5g"),
            if hit["station"] == bump_probe_station:
                hit["x"] -= self.target_bump[0]
                hit["px"] -= self.target_bump[1]
                print "move  ",
                self.target_hit = i
            elif hit["station"] in ignore_stations:
                print "ignore"
                continue
            else:
                print "      ",
            r_x = (hit["x"] - target_co[0])
            r_px = (hit["px"] - target_co[1])
            print format(r_x, "12.5g"),  format(r_px, "12.5g")
            x_score += r_x**2./n_hits
            px_score += r_px**2./n_hits
        return x_score, px_score


    def function(self, nvar, parameters, score, jacobian, err):
        self.iteration += 1
        x_scale = self.config.find_bump_parameters["position_tolerance"]
        px_scale = self.config.find_bump_parameters["momentum_tolerance"]
        print "Running minuit iteration", self.iteration, "target", self.target_bump
        fields = self.get_fields_from_minuit()
        print "    Fields", fields
        hit_list = self.track_one(fields)
        x_score, px_score = self.get_score(hit_list)
        print "Bump fields:"
        for i, bz in enumerate(fields):
            print "    bump", i+1, bz
        print "    x RMS residual ", format(x_score**0.5, "12.4g"), "scaled", format(x_score**0.5/x_scale, "12.4g")
        print "    px RMS residual", format(px_score**0.5, "12.4g"), "scaled", format(px_score**0.5/px_scale, "12.4g")
        score[0] = x_score/x_scale**2+px_score/px_scale**2
        self.score = score[0]
        print "score          ", format(score[0], "12.4g")
        print
        if self.iteration > self.max_iterations:
            raise RuntimeError("Hit maximum iteration")

    def setup_subs(self, fields):
        subs = self.config.substitution_list[0]
        if self.first_tracking:
            for key in sorted(subs.keys()):
                print utilities.sub_to_name(key), subs[key],
            self.first_tracking = False
        for key, value in self.overrides.iteritems():
            self.subs[key] = value
        for i, a_field in enumerate(fields):
            self.subs["__bump_field_"+str(i+1)+"__"] = a_field
        self.subs["__phi_foil_probe__"] = \
                             self.config.find_bump_parameters["foil_probe_phi"]
        xboa.common.substitute(self.config.tracking["lattice_file"],
                               self.tmp_dir+'SectorFFAGMagnet.tmp',
                               self.subs)

    def track_one(self, fields):
        ref_probes = self.config.find_bump_parameters["ref_probe_files"]
        #bump_probes = [self.config.find_bump_parameters["bump_probe_file"]]
        energy = self.config.find_bump_parameters["energy"]
        try:
            os.makedirs(self.tmp_dir)
        except OSError: # maybe the dir already exists
            pass
        self.setup_subs(fields)
        os.chdir(self.tmp_dir)
        find_closed_orbits.CONFIG = self.config
        ref_hit = find_closed_orbits.reference(energy)
        opal_exe = os.path.expandvars("${OPAL_EXE_PATH}/opal")
        tracking = OpalTracking('SectorFFAGMagnet.tmp', 'disttest.dat', ref_hit, ref_probes, opal_exe, "log")
        test_hit = ref_hit.deepcopy()
        closed_orbit = self.config.find_bump_parameters["closed_orbit"]
        test_hit["x"] = closed_orbit[0]
        test_hit["px"] = closed_orbit[1]
        # fix momentum
        test_hit["pz"] = (ref_hit["p"]**2-test_hit["px"]**2)**0.5
        print "Reference kinetic energy:", ref_hit["kinetic_energy"]
        print "Seed kinetic energy:     ", test_hit["kinetic_energy"]
        hit_list = tracking.track_one(test_hit)
        print "Station to probe mapping:\n   ",
        for i, fname in enumerate(tracking.get_names()):
            print "("+str(i)+",", fname+")",
        print

        self.tracking_result = [[hit["station"], hit["x"], hit["px"]] for hit in hit_list]
        return hit_list

def main(config):
    find_bump = FindBumpParameters(config)
    find_bump.find_bump_parameters()
    if __name__ == "__main__":
        raw_input()
