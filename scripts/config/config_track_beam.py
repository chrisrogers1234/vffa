import copy
import json
import glob
import math
import os
from . import config_arctan_triplet_baseline as config

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
            r0 = 0
        src_dir = os.path.join(os.getcwd(), 
                "output/arctan_baseline/bump_quest_v6/find_bump_r_"+str(r0)+"_theta_"+str(angle))
        self.run_control["output_dir"] = src_dir.replace("find", "track_")
        if float(r0) < 0:
            bump_fields = self.zero_bump_fields(config.get_baseline_substitution())
        else:
            bump_fields = self.softcoded_bump_fields(src_dir+"/find_bump_parameters_001.out")
        ring_tof = 1016.091866
        self.substitution_list = [config.get_baseline_substitution()] # used for plotting unperturbed CO
        self.substitution_list[0].update(bump_fields)
        self.run_control["find_closed_orbits_4d"] = True
        self.run_control["find_da"] = False
        self.run_control["find_bump_parameters"] = False
        self.run_control["track_bump"] = False
        self.run_control["track_beam"] = True
        self.find_closed_orbits["subs_overrides"]["__n_turns__"] = 1.11
        self.find_closed_orbits["final_subs_overrides"].update(bump_fields)
        self.find_closed_orbits["max_iterations"] = 0
        self.find_closed_orbits["do_minuit"] = True
        self.find_closed_orbits["minuit_iterations"] = 10      
        self.find_closed_orbits["us_cell"] = 0
        self.find_closed_orbits["ds_cell"] = 10

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

        max_amp = 0.01
        self.track_beam = {
            "run_dir":"tmp/track_beam/",
            "save_dir":"track_beam",
            "print_events":[0, 1, -1],
            "settings":[{
                "name":"forwards",
                "direction":"forwards",
                "probe_files":"RINGPROBE01.h5",          
                "beam":{
                    "type":"beam_gen",
                    "closed_orbit_file":"closed_orbits_cache",
                    "eigen_emittances":[[0, 0]]*2+[[max_amp, max_amp]],
                    "n_per_dimension":4,
                    "variables":["x","x'","y","y'"],
                    "amplitude_dist":"uniform", #"grid", # 
                    "phase_dist":"uniform", #"grid", # 
                    "max_amplitude_4d":max_amp, # amplitude_dist != grid
                    "energy":3.0,
                },
                "subs_overrides":{
                    "__n_turns__":26.1,
                    "__hdf5__":"True",
                    "__do_magnet_field_maps__":"False",
                    "__do_bump__":"True",
                    "__step_size__":0.01
                },
            },],
        }

        max_amp = 0.01
        self.track_beam = {
            "run_dir":"tmp/track_beam/",
            "save_dir":"track_beam",
            "print_events":[0, 1, -1],
            "settings":[{
                "name":"forwards",
                "direction":"forwards",
                "probe_files":"RINGPROBE01.h5",          
                "beam":{
                    "type":"beam_gen",
                    "closed_orbit_file":"closed_orbits_cache",
                    "eigen_emittances":[[0, 0]]*2+[[max_amp, max_amp]],
                    "n_per_dimension":4,
                    "variables":["x","x'","y","y'"],
                    "amplitude_dist":"uniform", #"grid", # 
                    "phase_dist":"uniform", #"grid", # 
                    "max_amplitude_4d":max_amp, # amplitude_dist != grid
                    "energy":3.0,
                },
                "subs_overrides":{
                    "__n_turns__":26.1,
                    "__hdf5__":"True",
                    "__do_magnet_field_maps__":"False",
                    "__do_bump__":"True",
                    "__step_size__":0.01
                },
            },],
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
        self.print_bump_fields(overrides)
        return overrides

    def zero_bump_fields(self, subs):
        for key in subs:
            if "bump" in key and "field" in key:
                subs[key] = 0.0
        subs["__do_bump__"] = True # Set to true to get foil probe...
        return subs



    def print_bump_fields(self, subs):
        print_str = ""
        print_dict = {}
        for key, value in subs.items():
            if "_bump_" in key and "_field__" in key:
                print_str += "    "+key.replace("__", "").replace("_", " ")+"   "+format(value, "8.5g")+"\n"
                print_dict[key] = subs[key]
        print(json.dumps(print_dict, indent=2))
        print(print_str)


if __name__ == "__main__":
    Config()
