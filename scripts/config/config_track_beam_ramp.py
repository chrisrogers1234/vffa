import copy
import json
import glob
import math
import os
from . import config_arctan_triplet_baseline as config

class Config(config.Config):
    def __init__(self):
        """
        Tracking to mimic "single turn injection" at low amplitude
        """
        super(Config, self).__init__()
        self.run_control["output_dir"] = os.getcwd()+"/output/arctan_baseline/single_turn_injection/tracking"
        ring_tof = 1149.185123
        field_files = glob.glob("output/arctan_baseline/bump_quest_v9/find_bump_r_*_theta_90/find_bump_parameters_001.out")
        substitutions = config.get_baseline_substitution(0.22, ring_tof)
        substitutions = self.ramp_fields(
            substitutions,
            field_files,
            [i for i, f in enumerate(field_files)],
            ring_tof,
            will_step = False
        )
        self.substitution_list = [substitutions] # used for plotting unperturbed CO

        self.run_control["find_closed_orbits_4d"] = False
        self.run_control["find_da"] = False
        self.run_control["find_bump_parameters"] = False
        self.run_control["track_bump"] = False
        self.run_control["track_beam"] = True
        self.find_closed_orbits["subs_overrides"]["__n_turns__"] = 0.11
        self.find_closed_orbits["subs_overrides"]["__do_bump__"] = False
        self.find_closed_orbits["final_subs_overrides"].update(substitutions)
        self.find_closed_orbits["max_iterations"] = 0
        self.find_closed_orbits["do_minuit"] = True
        self.find_closed_orbits["minuit_iterations"] = 10     
        self.find_closed_orbits["us_cell"] = 0
        self.find_closed_orbits["ds_cell"] = 1

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
        self.track_beam_dummy = {
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
                    "eigen_emittances":[[0, 0]]*3,#+[[max_amp, max_amp]],
                    "n_per_dimension":2,
                    "variables":["x","x'","y","y'"],
                    "amplitude_dist":"uniform", #"grid", # 
                    "phase_dist":"uniform", #"grid", # 
                    "max_amplitude_4d":max_amp, # amplitude_dist != grid
                    "energy":3.0,
                },
                "subs_overrides":{
                    "__n_turns__":45.1,
                    "__hdf5__":"True",
                    "__do_magnet_field_maps__":"False",
                    "__do_bump__":"True",
                    "__step_size__":0.01
                },
            },],
        }

        T0 = ring_tof
        self.track_beam = {
            "run_dir":"tmp/track_beam/",
            "save_dir":"track_beam_rf_on_2",
            "print_events":[i for i in range(1)],
            "variables":["x", "x'", "y", "y'", "t", "energy"],
            "settings":[{
                "name":"grid",
                "direction":"forwards",
                "probe_files":"RINGPROBE01.h5",          
                "beam":{
                    "type":"grid",
                    "energy":3.0,
                    "start":[4357.646683446333, 0.0, -116.7090485272821, 0.0, 0.0, 941.272],
                    "stop":[4357.646683446333, 0.0, -116.7090485272821, 0.0, T0, 941.272],
                    "nsteps":[1, 1, 1, 1, 4+1, 1],
                },
                "subs_overrides":{
                    "__n_turns__":100.1,
                    "__hdf5__":True,
                    "__do_magnet_field_maps__":False,
                    "__do_bump__":True,
                    "__do_rf__":True,
                    "__do_foil__":False,
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

    def sort_list(self, file_list):
        file_list_tmp = []
        for name in file_list:
            r = name.split("find_bump_r_")[1]
            r = r.split("_theta")[0]
            r = float(r)
            theta = name.split("_theta_")[1]
            theta = theta.split("/find_bump")[0]
            theta = float(theta)
            file_list_tmp.append([(r, theta), name])
        file_list_tmp = sorted(file_list_tmp)
        file_list_tmp = [f for f in reversed(file_list_tmp)]
        return file_list_tmp

    def ramp_fields(self, subs, file_list, turn_list, ring_tof, will_step):
        json_data_list = []
        for nominal_coordinates, file_name in self.sort_list(file_list):
            fin = open(file_name)
            best = None
            for line in fin.readlines():
                json_line = json.loads(line)
                if not best or json_line["score"] < best["score"]:
                    best = json_line
            print("Adding beam at", nominal_coordinates, "from file", file_name, "with score", best["score"], "Tracking", [t for t in best["tracking"] if t[0] == 4] )
            json_data_list.append(best)

        overrides = {}
        for file_index, json_line in enumerate(json_data_list):
            for key, value in json_line["subs"].items():
                if "_bump_" not in key or "_field_" not in key:
                    continue
                scale_key = key.replace("field", "scales")
                time_key = key.replace("field", "times")
                azimuth_key = key.replace("field", "phi")
                azimuth = subs[azimuth_key]/subs["__n_cells__"]
                if scale_key not in overrides:
                    overrides[key] = 1.0 #scale from 1.0
                    overrides[scale_key] = [value]
                    overrides[time_key] = [-1e9]
                if will_step:
                    overrides[scale_key] += [value, value]
                    overrides[time_key] += [(file_index-0.25)*ring_tof+1e-9, (file_index+0.75)*ring_tof]
                else:
                    time = (turn_list[file_index]+azimuth)*ring_tof #there is a phase shift owing to the azimuthal angle of the particle
                    overrides[scale_key] += [value]
                    overrides[time_key] += [time]
        for key in overrides.keys():
            if "bump" in key and "scales" in key:
                overrides[key] += [0.0, 0.0]
            elif "bump" in key and "times" in key:
                overrides[key] += [overrides[key][-1]+ring_tof, 1e9]
        subs.update(overrides)
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
