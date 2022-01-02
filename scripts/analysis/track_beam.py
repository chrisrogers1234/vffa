import os
import json
import copy
import math
import shutil
import glob
import time

import numpy
import xboa.hit

import plotting.plot_fields

from opal_tracking import OpalTracking
import utils.utilities
from utils.decoupled_transfer_matrix import DecoupledTransferMatrix
DecoupledTransferMatrix.det_tolerance = 1.


class BeamGen(object):
    def __init__(self, config):
        self.config = config

    def load_closed_orbit(self):
        file_name = os.path.join(self.config.run_control["output_dir"],
                                 self.config.track_beam["closed_orbit_file"])
        fin = open(file_name)
        self.closed_orbits  = json.loads(fin.read())

    def gen_beam(self):
        raise NotImplementedError("Not implemented")

    @classmethod
    def beam_setting(cls, beam, config):
        #if setting["beam"]["type"] == "last":
        #    station = setting["beam"]["station"]
        #    hit_list = [hit_list[station] for hit_list in self.last_tracking if station < len(hit_list)]
        if beam["type"] == "beam_gen":
            beam_gen = BeamShells(config, beam)
            hit_list = beam_gen.gen_beam()
        elif beam["type"] == "grid":
            beam_gen = BeamGrid(config, beam)
            hit_list = beam_gen.gen_beam()
        elif beam["type"] == "multibeam":
            beam_gen = MultiBeam(config, beam)
            hit_list = beam_gen.gen_beam()
        else:
            raise ValueError("Did not recognise beam of type "+str(setting["beam"]["type"]))
        return hit_list

class MultiBeam(object):
    def __init__(self, config, beam):
        self.config = copy.deepcopy(config)
        self.beam = beam
        self.beam_list = []

    def gen_beam(self):    
        hit_list = []
        for beam in self.beam["beam_list"]:
            hit_list += BeamGen.beam_setting(beam, self.config)
        return hit_list


class BeamShells(BeamGen):
    def __init__(self, config, beam):
        self.config = copy.deepcopy(config)
        self.config.track_beam = beam
        self.tm = None
        self.mean = None
        self.betagamma = 1.0

    def get_mean(self, co):
        hit = xboa.hit.Hit.new_from_dict(co["seed_hit"])
        self.mean = [hit[var] for var in self.config.track_beam["variables"]]
        self.betagamma = hit["p"]/hit["mass"]

    def get_tm(self, co):
        tm = co["tm"]
        tm = copy.deepcopy(co["tm"])
        for i, row in enumerate(tm):
            tm[i] = row[1:5]
        self.tm = DecoupledTransferMatrix(tm, True)

    def get_action(self, actions):
        dist = self.config.track_beam["amplitude_dist"]
        min_amp = 0.0
        max_amp = self.config.track_beam["max_amplitude_4d"]
        if dist == "grid":
            return actions[0], actions[1]
        au, av = -1, -1
        while au+av > max_amp or au+av < min_amp:
            if dist == "uniform":
                au = numpy.random.uniform(min_amp, actions[0])
                av = numpy.random.uniform(min_amp, actions[1])
            else:
                raise KeyError("Did not recognise amplitude_dist type "+str(dist))
        return au, av

    def get_phi(self, actions):
        n_phi_u = self.config.track_beam["n_per_dimension"]
        n_phi_v = self.config.track_beam["n_per_dimension"]
        if actions[0] == 0.0:
            n_phi_u = 1
        if actions[1] == 0.0:
            n_phi_v = 1
        phi_list = []
        dist = self.config.track_beam["phase_dist"]
        for ui in range(n_phi_u):
            for vi in range(n_phi_v):
                if dist == "grid":
                    phi_list.append([ui*2.*math.pi/n_phi_u-math.pi,
                                     vi*2.*math.pi/n_phi_v-math.pi])
                elif dist == "uniform":
                    phi_list.append([numpy.random.uniform(-math.pi, math.pi),
                                     numpy.random.uniform(-math.pi, math.pi)])
        return phi_list

    def get_aa_list(self, actions):
        aa_list = []
        phi_list = self.get_phi(actions)
        for phi_u, phi_v in phi_list:
            au, av = self.get_action(actions)
            aa = [phi_u, au, phi_v, av]
            aa_list.append(aa)
        return aa_list

    def get_psv_list(self, aa_list):
        psv_list = []
        for aa in aa_list:
            coupled = self.tm.action_angle_to_coupled(aa)
            psv = [coupled[i]+self.mean[i] for i in range(4)]
            psv_list.append(psv)
            print("aa:", aa, "coupled:", coupled, "psv:", psv)
        return psv_list

    def get_hit_list(self, psv_list):
        hit_list = []
        print("Tracking following particles:")
        for psv in psv_list:
            print(psv)
            mass = xboa.common.pdg_pid_to_mass[2212]
            energy = self.config.track_beam["energy"]+mass
            pz = (energy**2-mass**2)**0.5
            hit_dict = {"energy":energy, "mass":mass, "pid":2212, "pz":pz}
            for i, var in enumerate(self.config.track_beam["variables"]):
                hit_dict[var] = psv[i]
            hit = xboa.hit.Hit.new_from_dict(hit_dict, "pz")
            hit_list.append(hit)
        return hit_list

    def gen_beam(self):
        self.load_closed_orbit()
        hit_list = []
        for co in self.closed_orbits:
            self.get_mean(co)
            self.get_tm(co)
            eigen_emittance_list = self.config.track_beam["eigen_emittances"]
            for actions in eigen_emittance_list:
                aa_list = self.get_aa_list(actions) # action angle
                for aa in aa_list:
                    aa[1] /= self.betagamma
                    aa[3] /= self.betagamma
                psv_list = self.get_psv_list(aa_list) # phase space vector
                hit_list += self.get_hit_list(psv_list) # hit
        return hit_list

class BeamGrid(BeamGen):
    def __init__(self, config, beam):
        self.config = config
        self.beam = beam
        self.dim = len(beam["start"])
        self.reference = beam["reference"]
        self.start = beam["start"]
        self.stop = beam["stop"]
        self.nsteps = beam["nsteps"]
        self.step = [0.0 for i in range(self.dim)]
        self.hit_list = []
        for i in range(self.dim):
            if self.nsteps[i] == 1:
                continue
            self.step[i] = (beam["stop"][i]-beam["start"][i])/(beam["nsteps"][i]-1)
        
    def gen_grid(self):
        point = [0, 0, 0, 0, 0, 0]
        self.hit_list.append(copy.deepcopy(point))
        while True:
            for i in range(self.dim-1):
                if point[i] == self.nsteps[i]:
                    point[i] = 0
                    point[i+1] += 1
            if point[self.dim-1] == self.nsteps[self.dim-1]:
                break         
            self.hit_list.append(copy.deepcopy(point))
            point[0] += 1

    def gen_hits(self):
        print("Tracking following particles:")
        hit_list = []
        for point in self.hit_list:
            print("Adding", point)
            mass = xboa.common.pdg_pid_to_mass[2212]
            energy = self.beam["energy"]+mass
            pz = (energy**2-mass**2)**0.5
            hit_dict = {"energy":energy, "mass":mass, "pid":2212, "pz":pz}
            for i in range(self.dim):
                point[i] = self.start[i]+point[i]*self.step[i]
            for i, var in enumerate(self.config.track_beam["variables"]):
                hit_dict[var] = point[i]
            hit = xboa.hit.Hit.new_from_dict(hit_dict, "pz")
            hit_list.append(hit)
        if self.reference:
            ref_dict = {"energy":energy, "mass":mass, "pid":2212, "pz":pz}
            for i, var in enumerate(self.config.track_beam["variables"]):
                ref_dict[var] = self.reference[i]
            ref_hist = xboa.hit.Hit.new_from_dict(ref_dict, "pz")
            hit_list = [ref_hist]+hit_list
        self.hit_list = hit_list

    def gen_beam(self):
        self.gen_grid()
        self.gen_hits()
        return self.hit_list


class TrackBeam(object):
    def __init__(self, config):
        self.config = config
        self.tmp_dir = os.path.join(self.config.run_control["output_dir"],
                               self.config.track_beam["run_dir"])
        self.last_tracking = None

    def setup(self):
        #self.beam_gen = BeamShells(self.config)
        self.hit_list = None #self.beam_gen.gen_beam()

    def save_tracking(self, out_dir):
        src_dir = self.tmp_dir
        base_dir = os.path.join(self.config.run_control["output_dir"],
                                self.config.track_beam["save_dir"])
        target_dir = os.path.join(base_dir, out_dir)
        print("Saving to", target_dir)
        if not os.path.exists(base_dir):
            os.makedirs(base_dir)
        if os.path.exists(target_dir):
            shutil.rmtree(target_dir)
        shutil.copytree(src_dir, target_dir)
        time.sleep(1)

    def do_tracking(self):
        if len(self.config.substitution_list) > 1:
            raise RuntimeError("Didnt code subs list > 1")
        utils.utilities.clear_dir(self.tmp_dir)
        here = os.getcwd()
        os.chdir(self.tmp_dir)

        self.config.tracking["analysis_coordinate_system"] = "opal"
        for setting in self.config.track_beam["settings"]:
            self.track_setting(setting)
        print("done")
        os.chdir(here)

    def beam_setting(self, setting):
        if setting["beam"]["type"] == "last":
            station = setting["beam"]["station"]
            self.hit_list = [hit_list[station] for hit_list in self.last_tracking if station < len(hit_list)]
        else:
            self.hit_list = BeamGen.beam_setting(setting["beam"], self.config)

    def track_setting(self, setting): 
        print("Starting setting", setting["name"], "#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#")      
        self.tracking = utils.utilities.setup_tracking(self.config,
                                                  setting["probe_files"],
                                                  setting["beam"]["energy"])
        utils.utilities.do_lattice(self.config,
                                   self.config.substitution_list[0],
                                   setting["subs_overrides"])
        self.beam_setting(setting)
        if setting["direction"] == "backwards":
            for hit in self.hit_list:
                for var in "px", "py", "pz":
                    hit[var] *= -1
        elif setting["direction"] != "forwards":
            raise RuntimeError("Direction must be forwards or backwards")

        print("    ... tracking", len(self.hit_list), "hits", setting["direction"])
        if self.config.track_beam["print_events"] == "all":
            print_events = [i for i in range(len(self.hit_list))]
        else:
            print_events = self.config.track_beam["print_events"]
        for i in print_events:
            print("      Event", i, end="  ")
            for var in ["x", "y", "z", "px", "py", "pz", "t"]:
                print(var+":", self.hit_list[i][var], end=" ")
            print()

        self.last_tracking = self.tracking.track_many(self.hit_list)

        print("    ... found", [len(hit_list) for hit_list in self.last_tracking], "output hits")
        print("    from", self.tracking.get_name_dict())

        self.save_tracking(setting["name"])

        for filename in glob.glob("*.h5"):
            os.unlink(filename)
        print()

def main(config):
    tracker = TrackBeam(config)
    tracker.setup()
    tracker.do_tracking()

"""
Tracking 2 tracks
    ... forwards -  2 hits
x: 3778.0377273690738 y: -103.98094256513187 z: 0.0 px: 0.0 py: 0.0 pz: 75.09082484564964 
Saving to /mnt/db89a8c2-6f6b-4c5c-bbc4-5534fb22ff7a/cr67/work/2017-07-07_isis2/vertical_isis2/output/double_triplet_baseline/single_turn_injection/track_bump_parameters_x_10.0_y_0.0_mm_3/track_beam/forwards
    ... backwards -  2 hits
x: 3778.037727369074 y: -8.470329472543003e-19 z: -103.98094256513187 px: -2.5431913411137762e-17 py: -75.0908214268802 pz: 1.5894945881961101e-18 
"""

"""
None:
Tracking 2 tracks
    ... forwards -  2 hits
x: 3778.0377273690738 y: -103.98094256513187 z: 0.0 px: 0.0 py: 0.0 pz: 75.09082484564964 
Saving to /mnt/db89a8c2-6f6b-4c5c-bbc4-5534fb22ff7a/cr67/work/2017-07-07_isis2/vertical_isis2/output/double_triplet_baseline/single_turn_injection/track_bump_parameters_x_10.0_y_0.0_mm_3/track_beam/forwards
    ... backwards -  2 hits
x: -8.470329472543003e-19 y: -103.98094256513187 z: 3778.037727369074 px: -75.0908214268802 py: 1.5894945881961101e-18 pz: -2.5431913411137762e-17 
"""