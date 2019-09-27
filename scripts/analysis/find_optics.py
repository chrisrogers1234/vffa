"""
Script to find the tune; drives xboa DPhiTuneFinder (FFT was making large side
bands which looked non-physical)
"""

import glob
import json
import sys
import math
import os
import shutil
import ROOT
from opal_tracking import OpalTracking
import xboa.common as common
from xboa.hit import Hit

from utils import utilities

class Optics(object):
    def __init__(self, config):
        """
        Find the tune. 

        -probe_file_name: name of a PROBE file from OPAL output, or None. If file_name
                    is specified, use that file_name in order to calculate tunes
                    (and generate plots), otherwise generate a new one by
                    tracking.
        -closed_orbits_file_name: name of a file containing closed orbits,
                    generated by e.g. 
        """
        self.closed_orbits_cached = None # filled by _load_closed_orbits
        self.tmp_dir = None
        self.unique_id = 1
        self.lattice_src = config.tracking["lattice_file"]
        self.probe_filename = config.find_optics["probe_files"]
        co_file = os.path.join(config.run_control["output_dir"], config.find_closed_orbits["output_file"]+".out")
        self._load_closed_orbits(co_file)
        self.delta_x = config.find_optics["delta_x"]
        self.delta_y = config.find_optics["delta_y"]
        self.do_axis = config.find_optics["axis"]
        self.opal = config.tracking["opal_path"]
        self.step_size = config.tracking["step_size"]
        self.config = config

        self.lattice = "/SectorFFAGMagnet.tmp"
        self.beam_file = "/disttest.dat"
        self.log_file = "/log"
        self.output_filename = config.find_optics["output_file"]
        self.output_dir = config.run_control["output_dir"]

    def find_tune_dphi(self):
        """
        Algorithm is to just calculate the turn-by-turn phase advance; this is
        done by evolving turn-by-turn the track; calculating a matched ellipse
        by looking at tracking output; transforming the ellipse into a circle
        using LU decomposition; then calculating the angle advanced.
        """
        cwd = os.getcwd()
        fout = open(self.output_filename, "w")
        index = 0
        for i, closed_orbit in enumerate(self.closed_orbits_cached):
            if self.row != None and i not in self.row:
                continue
            if len(closed_orbit["hits"]) == 0:
                print("Error - no closed orbit")
                continue
            index += 1
            if index >= self.config.find_tune["root_batch"]:
                ROOT.gROOT.SetBatch(True)
            subs = closed_orbit["substitutions"]
            for item, key in self.config.find_tune["subs_overrides"].items():
                subs[item] = key

            print("Finding optics with", end=' ') 
            for key in sorted(subs.keys()):
                print(utilities.sub_to_name(key), subs[key], end=' ')
            print()
            tune_info = {"substitutions":subs}
            for axis1, axis2, delta1, delta2 in [("x", "px", self.delta_x, 0.),
                                                 ("y", "py", self.delta_y, 0.)]:
                if self.do_axis != None and axis1 != self.do_axis:
                    continue
                hit = Hit.new_from_dict(closed_orbit["hits"][0])
                self._temp_dir()
                common.substitute(
                    self.lattice_src, 
                    self.tmp_dir+self.lattice,
                    subs)
                tracking = OpalTracking(self.tmp_dir+self.lattice,
                                        self.tmp_dir+self.beam_file,
                                        self._reference(hit),
                                        self.probe_filename,
                                        self.opal,
                                        self.tmp_dir+self.log_file)
                tracking.clear_path = self.tmp_dir+"/*.loss"

            for key in sorted(tune_info.keys()):
                if "signal" not in key and "dphi" not in key:
                    print("   ", key, tune_info[key])
            print(json.dumps(tune_info), file=fout)
            fout.flush()
        os.chdir(cwd)

    def _temp_dir(self):
        """Make a temporary directory for tune calculation"""
        tmp_dir = os.path.join(self.output_dir, "tmp/find_optics/")
        try:
            os.makedirs(tmp_dir)
        except OSError:
            pass
        os.chdir(tmp_dir)
        self.tmp_dir = "./"

    def _load_closed_orbits(self, filename):
        """Load closed orbits from a json file"""
        fin = open(filename)
        closed_orbits = [json.loads(line) for line in fin.readlines()]
        self.closed_orbits_cached = closed_orbits
        
        #closed_orbits_energy = [orbit[0] for orbit in closed_orbits]
        #closed_orbits_x = [orbit[1:][0] for orbit in closed_orbits]
        #closed_orbits_dict = dict(zip(closed_orbits_energy, closed_orbits_x))

    def _reference(self, hit):
        """Generate a reference particle"""
        hit = hit.deepcopy()
        hit["x"] = 0.
        hit["px"] = 0.
        return hit










