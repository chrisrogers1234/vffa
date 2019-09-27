"""
Script to find the RF set up; drives find closed orbit
"""

import os
import sys
import copy
import json
import math
import glob

import numpy
import ROOT

import xboa.common
from xboa.hit import Hit
sys.path.insert(1, "scripts")
from opal_tracking import OpalTracking
from utils import utilities

class TrackBump(object):
    def __init__(self, config):
        self.config = config
        self.tmp_dir = os.path.join(self.config.run_control["output_dir"],
                               self.config.track_bump["run_dir"])
        lattice_name = self.config.tracking["lattice_file_out"].split(".")[0]
        for item in glob.glob(self.tmp_dir+"/"+lattice_name+"-trackOrbit*"):
            try:
                print("Clearing", item)
                os.unlink(item)
            except OSError: # maybe the dir already exists
                pass
        try:
            os.makedirs(self.tmp_dir)
        except OSError: # maybe the dir already exists
            pass
        self.subs = {}
        self.tracking_result = []
        self.bump_data = []
        self.output = []
        self.injection_orbit = None

    def track_bump(self):
        os.chdir(self.tmp_dir)
        for self.subs in self.config.substitution_list:
            self.load_bump_parameters()
            for bumps in self.bump_data:
                bump_list = self.get_bump_list(bumps)
                self.back_track(bump_list)
                self.fore_track(bump_list)
                #self.track_painting(bump_list)

    def is_equal_list(self, list_1, list_2, tolerance):
        if len(list_1) != len(list_2):
            return False
        for i, item_1 in enumerate(list_1):
            item_2 = list_2[i]
            if abs(item_1-item_2) > tolerance:
                return False
        return True


    def get_bump_list(self, bumps):
        """
        Get a list of bump magnet settings and set up injection orbits
        - bumps: the loaded bump magnet settings
        Set up the bump settings defined in config.track_bump["bump_list"].
        First compare the loaded bump finder data with the desired bumps to
        find the appropriate bump finder data. Then select the injection orbit
        corresponding to the config.track_bump["bump_probe_station"], which is
        the probe where we consider the foil to be placed.
        
        The return value is a list, each item in the list comprising a tuple
        having
            bump_list[0] dict of field name:field strength
            bump_list[1] integer corresponding to the number of turns for which
                         that field would be held (to simulate painting)
            bump_list[2] closed orbit ?? position of the closed orbit for that
                         field setting on the bump probe ??
        We also set self.injection_orbit, corresponding to the closed orbit
        for which the injected H- sit on top of a particular bump
        """
        bump_list = []
        bump_input = self.config.track_bump["bump_list"]
        if bump_input == None:
            find_bump_list = self.config.find_bump_parameters["bump"]
            foil_co = self.config.find_bump_parameters["foil_closed_orbit"]
            co = self.config.find_bump_parameters["closed_orbit"]
            bump_input = []
            for bump in find_bump_list:
                bump_input.append([])
                for i, bump_i in enumerate(bump):
                    bump_input[-1].append(bump_i+foil_co[i]-co[i])
        if self.config.track_bump["n_turns_list"] == None:
            n_turns_list = [1 for i in bump_input]
        source_bump = None
        for i, bump_target in enumerate(bump_input):
            for item in bumps["bumps"]:
                if self.is_equal_list(item["target_bump"], bump_target, 1e-9):
                    source_bump = item["target_bump"]
                    break
                print(item["target_bump"], bump_target)
            if source_bump == None:
                raise RuntimeError("Failed to find bump for "+str(bump_target))
                #tracking_index = item["target_hit"]
            closed_orbit = [0., 0.]
            for row in item["tracking"]:
                if row[0] == self.config.track_bump["bump_probe_station"]:
                    injection = self.config.track_bump["injection_orbit"]
                    if bump_input[injection] == bump_target:
                        self.injection_orbit = row[1:]
                        break
            for row in item["tracking"]:
                if row[0] == self.config.track_bump["bump_probe_station"]:
                    closed_orbit = row[1:]
                    break
            bump_list.append((item["bump_fields"], n_turns_list[i], closed_orbit))
            print("Got bump list fields:", bump_list[0], "n turns", n_turns_list[i], "co", closed_orbit)
        if self.injection_orbit == None:
            print("Failed to find injection orbit element", end=' ')
            print(self.config.track_bump["injection_orbit"], end=' ')
            print("from bump list of length", len(bump_list))
        else:
            print("Updated injection orbit", self.injection_orbit)
        return bump_list

    def load_bump_parameters(self):
        fname = os.path.join(self.config.run_control["output_dir"],
                             self.config.track_bump["input_file"])
        str_in = open(fname).read()
        self.bump_data = json.loads(str_in)

    def get_filename_root(self):
        fname = self.config.run_control["output_dir"]+"/"
        fname += self.config.track_bump["output_file"]
        return fname

    def save_state(self, suffix):
        state = {
            "target_bump":self.target_bump,
            "bump_fields":self.get_fields_from_minuit(),
            "tracking":self.tracking_result,
            "n_iterations":self.iteration,
        }
        self.output[-1]["bumps"].append(state)
        fname = self.get_filename_root()+"."+suffix
        fout = open(fname, "w")
        print(json.dumps(self.output, indent=2), file=fout)

    def do_substitutions(self, fields, phi_init, charge):
        overrides = self.config.track_bump["subs_overrides"]
        #print json.dumps(self.subs, indent=2)
        overrides["__beam_phi_init__"] = phi_init
        self.subs["__beam_charge__"] = charge
        overrides.update(fields)
        utilities.do_lattice(self.config, self.subs, overrides)

    def fore_track(self, bump_list):
        file_index = 0
        for bump_fields, n_turns, closed_orbit in bump_list:
            foil_phi = self.subs["__foil_probe_phi__"]
            foil_probe = self.config.track_bump["foil_probe_files"]
            energy = self.config.track_bump["energy"]
            self.do_substitutions(bump_fields, foil_phi, +1.0)
            tracking = utilities.setup_tracking(self.config, foil_probe, energy)
            test_hit = utilities.reference(self.config, energy,
                                           closed_orbit[0], closed_orbit[1],
                                           closed_orbit[2], closed_orbit[3])
            hit_list = tracking.track_one(test_hit)
            lattice_name = self.config.tracking["lattice_file_out"].split(".")[0]
            file_index += 1
            os.rename(lattice_name+"-trackOrbit.dat",
                      lattice_name+"-trackOrbit-fore-"+str(file_index)+".dat")

    def back_track(self, bump_list):
        file_index = 0
        for bump_fields, n_turns, closed_orbit in bump_list:
            foil_phi = self.subs["__foil_probe_phi__"]
            injected_beam = self.injection_orbit
            foil_probe = self.config.track_bump["foil_probe_files"]
            energy = self.config.track_bump["energy"]
            self.do_substitutions(bump_fields, foil_phi, -1.0)
            tracking = utilities.setup_tracking(self.config, foil_probe, energy)
            test_hit = utilities.reference(self.config, energy,
                                           injected_beam[0], injected_beam[1],
                                           injected_beam[2], injected_beam[3])
            test_hit["px"] *= -1
            test_hit["py"] *= -1
            test_hit["pz"] *= -1
            hit_list = tracking.track_one(test_hit)
            lattice_name = self.config.tracking["lattice_file_out"].split(".")[0]
            file_index += 1
            os.rename(lattice_name+"-trackOrbit.dat",
                      lattice_name+"-trackOrbit-back-"+str(file_index)+".dat")

    def track_painting(self, bump_list):
        foil_probe = self.config.track_bump["foil_probe_files"]
        ramp_probe = self.config.track_bump["ramp_probe_files"]
        foil_phi = self.config.track_bump["foil_probe_phi"]
        ramp_phi = self.config.track_bump["ramp_probe_phi"]
        injected_beam = self.injection_orbit
        energy = self.config.track_bump["energy"]

        self.do_substitutions({}, 0, +1.0)
        foil_tracking = utilities.setup_tracking(self.config, foil_probe, energy)
        ramp_tracking = utilities.setup_tracking(self.config, ramp_probe, energy)
        foil_tracking.verbose = False
        ramp_tracking.verbose = False
        test_hit = utilities.reference(self.config, energy,
                                        injected_beam[0], injected_beam[1],
                                        injected_beam[2], injected_beam[3])
        print("Seed kinetic energy:     ", test_hit["kinetic_energy"])
        # fix momentum
        hit_list_in = []
        log = open("track_bump.log", "w")
        for bump_fields, n_turns, closed_orbit in bump_list:
            # ramp the fields
            print("Ramping fields", bump_fields, file=log)
            if len(hit_list_in) > 0:
                # track from the ramp point to the foil point
                print("Tracking to foil", file=log)
                self.do_substitutions(bump_fields, ramp_phi, +1.0)
                hit_list_of_lists = foil_tracking.track_many(hit_list_in)
                hit_list_in = [hit_list[-1] for hit_list in hit_list_of_lists]
            self.do_substitutions(bump_fields, foil_phi, +1.0)
            for i in range(n_turns):
                hit_list_in.append(copy.deepcopy(test_hit))
                hit_list_in = [copy.deepcopy(test_hit), copy.deepcopy(test_hit)]+hit_list_in
                hit_list_of_lists = foil_tracking.track_many(hit_list_in)
                hit_list_of_lists = hit_list_of_lists[2:]
                hit_list_in = []
                for hit_list in hit_list_of_lists:
                    if hit_list[-1]["t"] > 100.:
                        hit_list_in.append(hit_list[-1])
                self.output.append(hit_list_in)
                print("  After turn", i+1, "fields", end=' ', file=log)
                for field in bump_fields:
                    print(field, end=' ')
                print("found", len(hit_list_of_lists), "particles", file=log)
                for item in hit_list_in:
                    print("   ", \
                          format(item["station"], "6"), \
                          format(item["t"], "12.6g"), \
                          format(item["x"], "12.6g"), \
                          format(item["px"], "12.6g"), \
                          format(item["kinetic_energy"], "12.6g"), file=log)
            # now track from the foil to the ramp point
            print("Tracking to ramp point", file=log)
            hit_list_of_lists = ramp_tracking.track_many(hit_list_in)
            hit_list_in = [hit_list[-1] for hit_list in hit_list_of_lists]
        return hit_list

def main(config):
    track_bump = TrackBump(config)
    track_bump.track_bump()
    if __name__ == "__main__":
        input()
