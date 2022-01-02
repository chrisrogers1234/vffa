import copy
import math
import os
import xboa.common
from . import config_double_triplet_baseline as config

class Config(config.Config):
    def __init__(self, energy = None):
        super(Config, self).__init__()
        if energy == None:
            energy = 3.0
        else:
            energy = float(energy)

        co_3MeV = [3778.037835263286, 0.0, -103.98039823958965, 0.0]
        subs = config.get_baseline_substitution()
        subs["__energy__"] = energy
        subs["__step_size__"] = 0.01
        mass = xboa.common.pdg_pid_to_mass[2212]
        p1 = (-mass**2+(mass+energy)**2)**0.5
        p0 = (-mass**2+(mass+3.0)**2)**0.5
        m_value = subs["__m_index__"]
        delta_pos = math.log(p1/p0)/m_value*1e3
        co = copy.deepcopy(co_3MeV)
        co[2] += delta_pos
        print("Calculating CO for p0:", p0, "p1:", p1, "m:", m_value, "delta:", delta_pos)
        print("Closed orbit at 3.0 MeV:", co_3MeV)
        print("             at", energy, "MeV:", co)
        self.substitution_list = [subs]
        self.find_closed_orbits["seed"] = [co]
        self.run_control["output_dir"] = os.path.join(os.getcwd(), "output/double_triplet_baseline/energy_scaling/"+str(energy)+"/")
        self.run_control["find_closed_orbits_4d"] = True
        self.run_control["find_da"] = True
        self.run_control["find_bump_parameters"] = False
        self.run_control["track_bump"] = False
        self.find_closed_orbits["final_subs_overrides"]["__n_turns__"] = 1.1
        self.find_closed_orbits["max_iterations"] = 0

        self.find_da = {
            "run_dir":"tmp/find_da/",
            "probe_files":"RINGPROBE01.h5",
            "subs_overrides":{"__n_turns__":101.1, "__do_magnet_field_maps__":"False", "__step_size__":0.01},
            "get_output_file":"get_da",
            "scan_output_file":"scan_da",
            "row_list":None,
            "scan_x_list":[],
            "scan_y_list":[],
            "x_seed":10.,
            "y_seed":10.,
            "min_delta":0.9,
            "max_delta":1000.,
            "required_n_hits":100,
            "dt_tolerance":0.5, # fraction of closed orbit dt
            "max_iterations":3,
            "decoupled":False,
        }
