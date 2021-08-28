import copy
import json
import glob
import math
import os
from . import config_double_triplet_baseline as config

class Config(config.Config):
    def __init__(self, angle = None, r0 = None):
        """
        Tracking beam to try to understand the cause of reduced DA. Aim is to
        plot tune and Maxwellianness vs u amplitude
        """
        super(Config, self).__init__()
        if angle == None:
            angle = 0
        if r0 == None:
            r0 = 10
        angle = float(angle)
        r0 = float(r0)
        src_dir = os.path.join(os.getcwd(),
                "output/double_triplet_baseline/single_turn_injection/bump_scan/find_bump_r_"+str(r0)+"_theta_"+str(angle))
        self.run_control["output_dir"] = src_dir.replace("find", "track")+"_scan"
        bump_fields = self.softcoded_bump_fields(src_dir+"/find_bump_parameters_001.out")
        ring_tof = 1016.091866
        eigen_emittances = [[0, 0], [0, 0]]
        for x in [0.005, 0.01, 0.02]:
            for i in range(6):
                au = i*x/5
                av = x-au
                eigen_emittances += [[au, av]] 
        eigen_emittances = [[0, 0]]*2+[[0.0, 0.002*i] for i in range(1, 11)]#, [0.02, 0.02]]
        delta = 0.01
        self.substitution_list = [self.zero_bump_fields(config.get_baseline_substitution())] # used for plotting unperturbed CO
        self.run_control["find_closed_orbits_4d"] = True
        self.run_control["find_da"] = False
        self.run_control["find_bump_parameters"] = False
        self.run_control["track_bump"] = False
        self.run_control["track_beam"] = True
        self.find_closed_orbits["final_subs_overrides"]["__n_turns__"] = 1.11
        self.find_closed_orbits["final_subs_overrides"]["__step_size__"] = 0.01
        self.find_closed_orbits["final_subs_overrides"]["__output_algorithm__"] = "RK4"
        self.find_closed_orbits["subs_overrides"]["__n_turns__"] = 1.11
        self.find_closed_orbits["subs_overrides"]["__step_size__"] = 0.01
        self.find_closed_orbits["subs_overrides"]["__output_algorithm__"] = "RK4"
        self.find_closed_orbits["final_subs_overrides"].update(bump_fields)
        self.find_closed_orbits["plotting_subs"].update(bump_fields)
        self.find_closed_orbits["subs_overrides"].update(bump_fields)
        self.find_closed_orbits["max_iterations"] = 1
        self.find_closed_orbits["do_minuit"] = False # just want to find the new tunes
        self.find_closed_orbits["us_cell"] = 0
        self.find_closed_orbits["ds_cell"] = 10
        self.find_closed_orbits["seed"] = [[3778.0377273690738, 0.0, -103.98094256513187, 0.0]]
        self.find_closed_orbits["use_py_tracking"] = False
        self.find_closed_orbits["deltas"] = [delta, delta*0.1, delta, delta*0.1]
        self.find_closed_orbits["subs_overrides"]["__do_bump__"] = False #############

        self.find_da = {
            "run_dir":"tmp/find_da/",
            "probe_files":"RINGPROBE01.h5",
            "subs_overrides":{"__n_turns__":101.1, "__do_magnet_field_maps__":"False", "__step_size__":0.01},
            "get_output_file":"get_da",
            "scan_output_file":"scan_da",
            "row_list":None,
            "scan_x_list":[],
            "scan_y_list":[],
            "x_seed":1.,
            "y_seed":1.,
            "min_delta":0.9,
            "max_delta":1000.,
            "required_n_hits":100,
            "dt_tolerance":0.5, # fraction of closed orbit dt
            "max_iterations":15,
            "decoupled":True,
        }

        self.track_beam = {
            "run_dir":"tmp/track_beam/",
            "save_dir":"track_beam",
            "print_events":[0],
            "settings":[{
                "name":"forwards_"+str(emit),
                "direction":"forwards",
                "probe_files":"FOILPROBE_1.h5",          
                "beam":{
                    "type":"beam_gen",
                    "closed_orbit_file":"closed_orbits_cache",
                    "eigen_emittances":[[emit*0.002, 0.0]],
                    "n_per_dimension":1,
                    "variables":["x","x'","y","y'"],
                    "amplitude_dist":"grid", # "uniform"
                    "phase_dist":"grid", # "uniform"
                    "max_amplitude_4d":0.01, # amplitude_dist != grid
                    "energy":3.0,
                },
                "subs_overrides":{
                    "__n_turns__":100.1,
                    "__hdf5__":"True",
                    "__do_magnet_field_maps__":"False",
                    "__do_bump__":"True",
                    "__step_size__":0.01
                },
            } for emit in range(11)]
        }

    def softcoded_bump_fields(self, file_name):
        fin = open(file_name)
        last = {"score":None}
        for line in fin.readlines():
            json_data = json.loads(line)
            if last["score"] == None or last["score"] > json_data["score"]:
                last = json_data
        overrides = last["bump_fields"]
        overrides["__do_bump__"] = True
        return overrides

    def zero_bump_fields(self, subs):
        for key in subs:
            if "bump" in key and "field" in key:
                subs[key] = 0.0
        subs["__do_bump__"] = True
        return subs



if __name__ == "__main__":
    Config()
