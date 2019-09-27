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
from utils import utilities

class FindBumpParameters(object):
    def __init__(self, config):
        self.config = config
        self.tmp_dir = os.path.join(self.config.run_control["output_dir"],
                               self.config.find_bump_parameters["run_dir"])
        self.x_err = self.config.find_bump_parameters["position_tolerance"]
        self.p_err = self.config.find_bump_parameters["momentum_tolerance"]
        self.first_tracking = True
        self.field_names = []
        self.subs = {}
        self.iteration = 0
        self.target_bump = [0., 0., 0., 0.]
        self.tracking_result = []
        self.output = []
        self.score = None
        self.target_hit = None
        self.overrides = {}

    def setup_minuit(self):
        self.iteration = 0
        self.field_names = [name for name in self.field_names_gen()]
        self.minuit = ROOT.TMinuit(len(self.field_names))
        min_b = self.config.find_bump_parameters["magnet_min_field"]
        max_b = self.config.find_bump_parameters["magnet_max_field"]
        fixed = self.config.find_bump_parameters["fix_bumps"]
        seed = self.config.find_bump_parameters["seed_field"]
        errors = self.config.find_bump_parameters["seed_errors"]
        print(len(self.field_names), len(seed), len(errors))
        for i, name in enumerate(self.field_names):
            self.minuit.DefineParameter(i, name, seed[i], abs(errors[i]), min_b, max_b)
            if name in fixed:
                self.minuit.FixParameter(i)
                print("Fixed bump", name)

        self.minuit.SetFCN(self.function)

    def field_names_gen(self):
        n_h_bumps = self.config.find_bump_parameters["n_h_bumps"]
        n_v_bumps = self.config.find_bump_parameters["n_v_bumps"]
        for i in range(n_h_bumps):
            yield "__h_bump_"+str(i+1)+"_field__"
        for i in range(n_v_bumps):
            yield "__v_bump_"+str(i+1)+"_field__"

    def find_bump_parameters(self):
        for subs in self.config.substitution_list:
            self.output.append({"subs":subs, "bumps":[]})
            for bump in self.config.find_bump_parameters["bump"]:
                for i, bump_i in enumerate(bump):
                    foil_co = self.config.find_bump_parameters["foil_closed_orbit"]
                    co = self.config.find_bump_parameters["closed_orbit"]
                    self.target_bump[i] = bump_i+foil_co[i]-co[i]
                self.setup_minuit()
                self.subs = subs
                self.overrides = self.config.find_bump_parameters["subs_overrides"]
                self.max_iterations = self.config.find_bump_parameters["max_iterations"]
                try:
                    self.minuit.Command("SIMPLEX "+str(self.max_iterations)+" 1.")
                except Exception:
                    sys.excepthook(*sys.exc_info())
                    print("Minuit failed")
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
            "target_bump":copy.deepcopy(self.target_bump),
            "bump_fields":self.get_fields_from_minuit(),
            "tracking":copy.deepcopy(self.tracking_result),
            "n_iterations":self.iteration,
            "target_hit":copy.deepcopy(self.target_hit),
            "score":self.score,
        }
        if append:
            self.output[-1]["bumps"].append(state)
        fname = self.get_filename_root()+"."+suffix
        fout = open(fname, "w")
        print(json.dumps(self.output, indent=2), file=fout)
 
    def get_fields_from_minuit(self):
        fields = {}
        for i, name in enumerate(self.field_names):
            var = ROOT.Double()
            err = ROOT.Double()
            self.minuit.GetParameter(i, var, err)
            fields[name] = var
        return fields

    def get_score(self, hit_list):
        ignore_stations = self.config.find_bump_parameters["ignore_stations"]
        target_co = self.config.find_bump_parameters["closed_orbit"]
        bump_probe_station = self.config.find_bump_parameters["bump_probe_station"]
        min_n_hits = self.config.find_bump_parameters["min_n_hits"]
        n_hits = len(hit_list)-len(ignore_stations)
        x_score, px_score, y_score, py_score = 0., 0., 0., 0.
        for i, hit in enumerate(hit_list):
            print(format(hit["station"], "6"), end=' ')
            print(format(hit["x"], "8.5g"), format(hit["px"], "8.5g"), end=' ')
            print(format(hit["y"], "8.5g"), format(hit["py"], "8.5g"), end=' ')
            print("|", format(hit["t"], "12.5g"), end=' ')
            if hit["station"] == bump_probe_station:
                hit["x"] -= self.target_bump[0]
                hit["px"] -= self.target_bump[1]
                hit["y"] -= self.target_bump[2]
                hit["py"] -= self.target_bump[3]
                print("move  ", end=' ')
                self.target_hit = i
            elif hit["station"] in ignore_stations:
                print("ignore", "|")
                continue
            else:
                print("      ", end=' ')
            r_x = (hit["x"] - target_co[0])
            r_px = (hit["px"] - target_co[1])
            r_y = (hit["y"] - target_co[2])
            r_py = (hit["py"] - target_co[3])
            print("|", format(r_x, "8.3g"),  format(r_px, "8.3g"), end=' ')
            print(format(r_y, "8.3g"),  format(r_py, "8.3g"))
            x_score += r_x**2./n_hits
            px_score += r_px**2./n_hits
            y_score += r_y**2./n_hits
            py_score += r_py**2./n_hits
        if n_hits <= min_n_hits:
            print("Insufficient hits - penalising by factor 1e9!")
            x_score = x_score*1e9+1e9
            px_score = px_score*1e9+1e9
            y_score = y_score*1e9+1e9
            py_score = py_score*1e9+1e9
        return x_score, px_score, y_score, py_score


    def function(self, nvar, parameters, score, jacobian, err):
        self.iteration += 1
        x_scale = self.config.find_bump_parameters["position_tolerance"]
        px_scale = self.config.find_bump_parameters["momentum_tolerance"]
        z_scale = self.config.find_bump_parameters["position_tolerance"]
        pz_scale = self.config.find_bump_parameters["momentum_tolerance"]
        print("Running minuit iteration", self.iteration, "target", self.target_bump)
        fields = self.get_fields_from_minuit()
        hit_list = self.track_one(fields)
        x_score, px_score, z_score, pz_score = self.get_score(hit_list)
        print("Bump fields:")
        for key in sorted(fields.keys()):
            value = fields[key]
            print(key.replace("_", " "), value)
        print("    x RMS residual ", format(x_score**0.5, "12.4g"), "scaled", format(x_score**0.5/x_scale, "12.4g"))
        print("    px RMS residual", format(px_score**0.5, "12.4g"), "scaled", format(px_score**0.5/px_scale, "12.4g"))
        print("    z RMS residual ", format(z_score**0.5, "12.4g"), "scaled", format(z_score**0.5/z_scale, "12.4g"))
        print("    pz RMS residual", format(pz_score**0.5, "12.4g"), "scaled", format(pz_score**0.5/pz_scale, "12.4g"))
        score[0] = x_score/x_scale**2+px_score/px_scale**2+z_score/z_scale**2+pz_score/pz_scale**2
        self.score = score[0]
        print("score          ", format(score[0], "12.4g"))
        print()
        if self.iteration > self.max_iterations:
            raise RuntimeError("Hit maximum iteration")

    def setup_subs(self, fields):
        subs = self.config.substitution_list[0]
        if self.first_tracking:
            for key in sorted(subs.keys()):
                print(utilities.sub_to_name(key), subs[key], end=' ')
            self.first_tracking = False
        self.overrides.update(fields)
        print("OVERRIDES", self.overrides)
        utilities.do_lattice(self.config, self.subs, self.overrides)

    def track_one(self, fields):
        utilities.clear_dir(self.tmp_dir)
        os.chdir(self.tmp_dir)
        ref_probes = self.config.find_bump_parameters["ref_probe_files"]
        #bump_probes = [self.config.find_bump_parameters["bump_probe_file"]]
        ref_energy = self.config.find_bump_parameters["energy"]
        self.setup_subs(fields)
        tracking = utilities.setup_tracking(self.config, ref_probes, ref_energy)
        test_hit = tracking.ref.deepcopy()
        closed_orbit = self.config.find_bump_parameters["closed_orbit"]
        test_hit["x"] = closed_orbit[0]
        test_hit["px"] = closed_orbit[1]
        test_hit["y"] = closed_orbit[2]
        test_hit["py"] = closed_orbit[3]
        # fix momentum
        test_hit["pz"] = (tracking.ref["p"]**2-test_hit["px"]**2)**0.5
        print("Reference kinetic energy:", tracking.ref["kinetic_energy"])
        print("Seed kinetic energy:     ", test_hit["kinetic_energy"])
        hit_list = tracking.track_one(test_hit)
        print("Station to probe mapping:\n   ", end=' ')
        for i, fname in enumerate(tracking.get_names()):
            print("("+str(i)+",", fname+")", end=' ')
        print()

        self.tracking_result = [[hit["station"], hit["x"], hit["px"], hit["y"], hit["py"]] for hit in hit_list]
        return hit_list

def main(config):
    find_bump = FindBumpParameters(config)
    find_bump.find_bump_parameters()
    if __name__ == "__main__":
        input()
