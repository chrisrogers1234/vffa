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
import platypus

import xboa.common
from xboa.hit import Hit
sys.path.insert(1, "scripts")
from opal_tracking import OpalTracking
from utils import utilities
from utils.decoupled_transfer_matrix import DecoupledTransferMatrix

class FindBumpParameters(object):
    def __init__(self, config):
        self.config = config
        self.tmp_dir = os.path.join(self.config.run_control["output_dir"],
                               self.config.find_bump_parameters["run_dir"])
        self.x_err = self.config.find_bump_parameters["position_tolerance"]
        self.p_err = self.config.find_bump_parameters["momentum_tolerance"]
        self.field_names = [name for name in self.field_names_gen()]
        self.fixed = {}
        self.first_tracking = True
        self.subs = {}
        self.iteration = 0
        self.target_bump = [0., 0., 0., 0.]
        self.target_amp = [0., 0.]
        self.tracking_result = []
        self.output = []
        self.score = None
        self.target_hit = None
        self.overrides = {}
        self.bump_fields = {}
        self.store_index = 1
        DecoupledTransferMatrix.det_tolerance = 1.

    def setup_minuit(self):
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
        self.minuit.SetFCN(self.minuit_function)

    def setup_amplitude(self, substitution_index):
        tm_source = self.config.run_control["output_dir"]+"/"+\
                    self.config.find_bump_parameters["tm_source"]
        co_source = open(tm_source)
        co_source = json.loads(co_source.readline())
        tm = co_source[substitution_index]["tm"]
        for i, row in enumerate(tm):
            tm[i] = row[1:5]
        self.tm = DecoupledTransferMatrix(tm)

    def field_names_gen(self):
        n_h_bumps = self.config.find_bump_parameters["n_h_bumps"]
        n_v_bumps = self.config.find_bump_parameters["n_v_bumps"]
        for i in range(n_h_bumps):
            yield "__h_bump_"+str(i+1)+"_field__"
        for i in range(n_v_bumps):
            yield "__v_bump_"+str(i+1)+"_field__"

    def store_data(self):
        outname = self.get_filename_root()+"_"+str(self.store_index)+".out"
        tmpname = self.get_filename_root()+".tmp"
        print("Moving data from", tmpname, "to", outname)
        os.rename(tmpname, outname)
        self.store_index += 1

    def find_bump_parameters(self):
        try:
            os.rename(
                self.get_filename_root()+".tmp",
                self.get_filename_root()+".old"
            )
        except OSError:
            pass
        for i, subs in enumerate(self.config.substitution_list):
            self.setup_amplitude(i)
            self.output.append({"subs":subs, "bumps":[]})
            for bump in self.config.find_bump_parameters["bump"]:
                self.do_one_optimisation(bump, None, subs)
            for amp in self.config.find_bump_parameters["amp"]:
                self.do_one_optimisation(None, amp, subs)

    def run_minuit(self, algorithm):
        self.setup_minuit()
        try:
            self.minuit.Command("SIMPLEX "+str(self.max_iterations)+" 1.")
        except StopIteration:
            sys.excepthook(*sys.exc_info())
            print("Minuit failed")

    def run_platypus(self, algorithm):
        min_b = self.config.find_bump_parameters["magnet_min_field"]
        max_b = self.config.find_bump_parameters["magnet_max_field"]
        fixed = self.config.find_bump_parameters["fix_bumps"]
        seed = copy.deepcopy(self.config.find_bump_parameters["seed_field"])
        sigma = self.config.find_bump_parameters["seed_errors"][0]
        self.field_names = [name for name in self.field_names_gen()]
        max_iterations = self.config.find_bump_parameters["max_iterations"]
        n_target_fields = len(self.config.find_bump_parameters["target_fields"])


        for name in fixed:
            index = self.field_names.index(name)
            self.fixed[name] = seed[index]
            print("Fixing bump", name, "to", self.fixed[name])
            del self.field_names[index]
            del seed[index]
        n_parameters = len(self.field_names)
        n_variables = 6+n_target_fields
        problem = platypus.Problem(n_parameters, n_variables)
        problem.types[:] = platypus.Real(min_b, max_b)
        problem.function = self.moo_function
        algorithm = platypus.CMAES(problem)
        algorithm.initial_search_point = seed
        algorithm.sigma = sigma # note this is a single value
        try:
            algorithm.run(max_iterations)
        except Exception:
            sys.excepthook(*sys.exc_info())
            print("Platypus failed")

    def do_one_optimisation(self, bump, amp, subs):
        print("Doing optimisation with bump", bump, "amp", amp)
        if bump != None:
            for i, bump_i in enumerate(bump):
                foil_co = self.config.find_bump_parameters["foil_closed_orbit"]
                co = self.config.find_bump_parameters["closed_orbit"]
                self.target_bump[i] = bump_i+foil_co[i]-co[i]
        else:
            self.target_bump = None
        self.target_amp = amp
        self.subs = subs
        self.overrides = self.config.find_bump_parameters["subs_overrides"]
        self.max_iterations = self.config.find_bump_parameters["max_iterations"]
        self.iteration = 0
        algorithm = self.config.find_bump_parameters["algorithm"]
        if algorithm == "simplex":
            self.run_minuit(algorithm)
        elif algorithm == "nsga2":
            self.run_platypus(algorithm)
        else:
            raise RuntimeError("Did not recognise algorithm "+str(algorithm))
        print("Finished optimisation with bump", bump, "amp", amp)
        self.overrides = self.config.find_bump_parameters["final_subs_overrides"]
        print("Doing final tracking")
        self.track_one(self.bump_fields)
        print("Done final tracking")
        self.store_data()
        print("End of optimisation loop\n\n")

    def get_filename_root(self):
        fname = self.config.run_control["output_dir"]+"/"
        fname += self.config.find_bump_parameters["output_file"]
        return fname

    def save_state(self, suffix, append):
        state = {
            "target_bump":copy.deepcopy(self.target_bump),
            "target_amp":copy.deepcopy(self.target_amp),
            "bump_probe_station":self.config.find_bump_parameters["bump_probe_station"],
            "bump_fields":self.bump_fields,
            "tracking":copy.deepcopy(self.tracking_result),
            "n_iterations":self.iteration,
            "target_hit":copy.deepcopy(self.target_hit),
            "score":self.score,
            "subs":self.overrides,
        }
        fname = self.get_filename_root()+"."+suffix
        fout = open(fname, "a")
        print(json.dumps(state), file=fout)
 
    def get_fields_from_minuit(self):
        fields = {}
        for i, name in enumerate(self.field_names):
            var = ROOT.Double()
            err = ROOT.Double()
            self.minuit.GetParameter(i, var, err)
            fields[name] = var
        return fields
 
    def get_fields_from_parameters(self, parameters_list):
        fields = copy.deepcopy(self.fixed)
        for i, name in enumerate(self.field_names):
            par = parameters_list[i]
            fields[name] = par
        self.bump_fields = fields
        return fields

    def get_amplitude(self, hit):
        """
        Oh what a beautiful morning
        I have put in place the machinery, but I dont *quite* know what I am
        trying to optimise for... amplitude of the H- injection relative to the
        proton closed orbit; where we move the proton closed orbit and leave the
        H- injecting to a fixed point; TM at the foil needs to be evaluated.
        
        So I think I need to do:
            1. calculate 50 mm offset bumper settings
            2. calculate TM using the bumper settings (and e.g. CO finder)
            3. scan amplitude settings relative to bumper settings
        """
        return 0., 0.

    def get_n_hits(self, hit_list, ignore_stations):
        target_n_hits = self.config.find_bump_parameters["target_n_hits"]
        hit_list = [hit for hit in hit_list if hit["station"] not in ignore_stations]
        n_hits = len(hit_list)
        penalty_factor = self.config.find_bump_parameters["penalty_factor"]
        denominator = n_hits
        if denominator == 0:
            denominator = 1
        penalty = penalty_factor**(target_n_hits-n_hits)
        if penalty < 1:
            penalty = 1.
        print("Found ", n_hits, "usable hits.",
              "Target", target_n_hits, "hits.",
              "Penalty", penalty, "will be applied to all scores")
        return n_hits, denominator, penalty

    def get_score(self, hit_list):
        ignore_stations = self.config.find_bump_parameters["ignore_stations"]
        target_co = self.config.find_bump_parameters["closed_orbit"]
        bump_probe_station = self.config.find_bump_parameters["bump_probe_station"]
        amp_probe_station = self.config.find_bump_parameters["amp_probe_station"]
        x_score, px_score, y_score, py_score, au_score, av_score = 0., 0., 0., 0., 0., 0.
        n_hits, denominator, penalty = self.get_n_hits(hit_list, ignore_stations)
        for i, hit in enumerate(hit_list):
            print(format(hit["station"], "6"), end=' ')
            print(format(hit["x"], "12.9g"), format(hit["px"], "8.5g"), end=' ')
            print(format(hit["y"], "8.5g"), format(hit["py"], "8.5g"), end=' ')
            print("|", format(hit["t"], "12.5g"), end=' ')
            au, av = 0., 0.
            if hit["station"] == bump_probe_station and self.target_bump != None:
                hit["x"] -= self.target_bump[0]
                hit["px"] -= self.target_bump[1]
                hit["y"] -= self.target_bump[2]
                hit["py"] -= self.target_bump[3]
                print("move  ", end=' ')
                self.target_hit = i
            elif hit["station"] in ignore_stations:
                print("ignore", "|")
                continue
            elif hit["station"] == amp_probe_station and self.target_amp != None:
                au, av = self.get_amplitude(hit)
                hit["x"] = target_co[0]
                hit["px"] = target_co[1]
                hit["y"] = target_co[2]
                hit["py"] = target_co[3]
                print("amp  ", end=' ')
            else:
                print("      ", end=' ')
            r_x = (hit["x"] - target_co[0])
            r_px = (hit["px"] - target_co[1])
            r_y = (hit["y"] - target_co[2])
            r_py = (hit["py"] - target_co[3])
            print("|", format(r_x, "8.3g"),  format(r_px, "8.3g"), end=' ')
            print(format(r_y, "8.3g"),  format(r_py, "8.3g"))
            x_score += r_x**2./denominator
            px_score += r_px**2./denominator
            y_score += r_y**2./denominator
            py_score += r_py**2./denominator
            au_score += au**2/denominator
            av_score += av**2/denominator
        x_score = x_score*penalty
        px_score = px_score*penalty
        y_score = y_score*penalty
        py_score = py_score*penalty
        return x_score, px_score, y_score, py_score, au_score, av_score

    def minuit_function(self, nvar, parameters, score, jacobian, err):
        fields = self.get_fields_from_minuit()
        # a bit stupid, to get the interface right we convert from minuit to 
        # dict to list to dict
        parameters = [fields[name] for name in self.field_names]
        score_list = self.moo_function(parameters)
        score[0] =  self.score

    def moo_function(self, parameters):
        self.iteration += 1
        if self.iteration > self.max_iterations:
            raise StopIteration("Hit maximum iteration")
        pos_scale = self.config.find_bump_parameters["position_tolerance"]
        mom_scale = self.config.find_bump_parameters["momentum_tolerance"]
        a_scale = self.config.find_bump_parameters["amplitude_tolerance"]
        b_scale = self.config.find_bump_parameters["field_tolerance"]
        target_fields = self.config.find_bump_parameters["target_fields"]

        print("Running iteration", self.iteration, "target", self.target_bump)
        fields = self.get_fields_from_parameters(parameters)
        hit_list = self.track_one(fields)
        x_score, px_score, z_score, pz_score, au_score, av_score = self.get_score(hit_list)
        b_score = [(fields[key] - target_fields[key])**2/b_scale**2 for key in target_fields ]
        print("Bump fields:")
        print("   ", [fields[key] for key in sorted(fields)])
        for key in sorted(fields.keys()):
            value = fields[key]
            print(key.replace("_", " "), value)
        score_list = [x_score/pos_scale**2,
               px_score/mom_scale**2,
               z_score/pos_scale**2,
               pz_score/mom_scale**2,
               au_score/a_scale**2,
               av_score/a_scale**2]+b_score
        self.score = sum(score_list)

        print("Residuals are RMS*(penalty due to missing hits)")
        print("    x RMS residual ", format(x_score**0.5, "12.4g"))
        print("    px RMS residual", format(px_score**0.5, "12.4g"))
        print("    z RMS residual ", format(z_score**0.5, "12.4g"))
        print("    pz RMS residual", format(pz_score**0.5, "12.4g"))
        print("    au RMS residual", format(au_score**0.5, "12.4g"))
        print("    av RMS residual", format(av_score**0.5, "12.4g"))
        print("    b score        ", format(sum(b_score)**0.5*b_scale, "12.4g"))
        print("score          ", format(self.score, "12.4g"))
        print()
        self.save_state("tmp", True)
        return score_list

    def setup_subs(self, fields):
        subs = self.config.substitution_list[0]
        if self.first_tracking:
            for key in sorted(subs.keys()):
                print(utilities.sub_to_name(key), subs[key], end=' ')
            self.first_tracking = False
        self.overrides.update(fields)
        print("OVERRIDES", self.overrides)
        utilities.do_lattice(self.config, self.subs, self.overrides)

    def cuts(self, hit_list):
        min_delta = self.config.find_bump_parameters["min_time_delta"]
        index = 0
        while index+1 < len(hit_list):
            if hit_list[index+1]["t"] - hit_list[index]["t"] < min_delta:
                del hit_list[index+1]
            elif hit_list[index+1]["station"] == hit_list[index]["station"]:
                del hit_list[index+1]
            else:
                index += 1
        return hit_list

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
        hit_list = tracking.track_many([test_hit])[1]
        print("Station to probe mapping:\n   ", end=' ')
        for i, fname in enumerate(tracking.get_names()):
            print("("+str(i)+",", fname+")", end=' ')
        print()
        hit_list = self.cuts(hit_list)
        self.tracking_result = [[hit["station"], hit["x"], hit["px"], hit["y"], hit["py"]] for hit in hit_list]
        return hit_list

def main(config):
    find_bump = FindBumpParameters(config)
    find_bump.find_bump_parameters()
    if __name__ == "__main__":
        input()
