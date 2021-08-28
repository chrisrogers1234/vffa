import copy
import json
import glob
import math
import os
from . import config_double_triplet_baseline as config

class Config(config.Config):
    def __init__(self):

        ring_tof = 1016.091866
        ramp_turns = 1
        bump_fields = {}
        if ramp_turns != None:
            end_time = ramp_turns*ring_tof
            bump_fields = self.get_bump_subs(0, end_time, ring_tof)

        eigen_emittances = [[0, 0], [0, 0]]
        for x in [1e-3]:#, 2e-3, 4e-3, 8e-3]:
            eigen_emittances += [[0, x], [x, 0], [x, x]] 
        eigen_emittances = [[0,0], [0,0], [0,0], [0, 1e-4]]
        super(Config, self).__init__()
        delta = 0.01
        self.substitution_list = [config.get_baseline_substitution()]
        self.run_control["find_closed_orbits_4d"] = True
        self.run_control["find_da"] = False
        self.run_control["find_bump_parameters"] = False
        self.run_control["track_bump"] = False
        self.run_control["track_beam"] = False
        self.find_closed_orbits["final_subs_overrides"]["__n_turns__"] = 1.51
        self.find_closed_orbits["final_subs_overrides"]["__step_size__"] = 0.001
        self.find_closed_orbits["final_subs_overrides"]["__output_algorithm__"] = "RK4"
        self.find_closed_orbits["final_subs_overrides"].update(bump_fields)
        self.find_closed_orbits["plotting_subs"].update(bump_fields)
        self.find_closed_orbits["subs_overrides"]["__n_turns__"] = 0.61
        self.find_closed_orbits["subs_overrides"]["__step_size__"] = 0.001
        self.find_closed_orbits["subs_overrides"]["__output_algorithm__"] = "RK4"
        self.find_closed_orbits["max_iterations"] = 0
        self.find_closed_orbits["us_cell"] = 4
        self.find_closed_orbits["ds_cell"] = 5
        self.find_closed_orbits["seed"] = [[3792.92652, +2.21609692/75, -63.9803651, -0.938736658/75]]
        self.find_closed_orbits["use_py_tracking"] = False
        self.find_closed_orbits["deltas"] = [delta, delta*0.1, delta, delta*0.1]
        self.find_closed_orbits["do_minuit"] = False

        self.track_beam = {
            "closed_orbit_file":"closed_orbits_cache",
            "eigen_emittances":eigen_emittances,# [1, 0], [1, 1]],
            "n_per_dimension":10,
            "variables":["x","x'","y","y'"],
            "energy":3.0,
            "probe_files":"RINGPROBE*.h5",
            "subs_overrides":{
                "__n_turns__":1.1,
                "__hdf5__":"True",
                "__do_magnet_field_maps__":"False",
            },
            "run_dir":"tmp/track_beam/",

        }
        self.track_beam["subs_overrides"].update(bump_fields)

        self.run_control["output_dir"] = os.path.join(os.getcwd(), "output/double_triplet_baseline/kicker_injection/track_bump_1")

    def get_bump_subs(self, start_time, end_time, ring_tof):
        subs = {"__do_bump__":True}
        bump_fields = self.hardcoded_bump_fields()
        if len(bump_fields):
            time_step = (end_time-start_time)/len(bump_fields)
        else:
            return {}
        for bump_field in reversed(bump_fields):
            # list of field values __h_bump_1_field__
            # list of scalings     __h_bump_1_scales__ 
            # list of times        __h_bump_1_times__
            for key in bump_field.keys():
                scale_key = key.replace("field", "scales")
                if scale_key not in subs:
                    subs[scale_key] = [bump_field[key]] # we duplicate the first value to model the flat top
                subs[scale_key].append(bump_field[key])
        times = [-1e9]+[start_time+time_step*i for i in range(len(bump_fields)+1)]+[1e9]
        time_delays = self.get_bump_time_delays(ring_tof)
        for key in bump_field.keys():
            scale_key = key.replace("field", "scales")
            subs[scale_key].append(bump_field[key]) # we also triplicate the last value to model the flat bottom
            subs[scale_key].append(bump_field[key])
            delay_key = key.replace("field", "delay")
            my_times = [t+time_delays[delay_key] for t in times]
            time_key = key.replace("field", "times")
            subs[time_key] = my_times
            subs[key] = 1.0
            print(len(subs[time_key]), len(subs[scale_key]))
        for key in sorted(subs.keys()):
            print(key, subs[key])
        return subs

    def get_bump_time_delays(self, ring_tof):
        subs = config.get_baseline_substitution()
        time_delays = {}
        for key in subs.keys():
            if "bump" not in key or "phi" not in key:
                continue
            delay_key = key.replace("phi", "delay")
            time_delays[delay_key] = subs[key]*ring_tof/10.0
        return time_delays


    def get_bump_fields(self):
        file_name_glob = "output/double_triplet_baseline/single_turn_injection/find_bump_parameters_4/find_bump_parameters*.out"
        file_name_list = sorted(glob.glob(file_name_glob))
        file_name_list = [fname for i, fname in enumerate(file_name_list)]

        field_list = []
        for file_name in file_name_list:
            fin = open(file_name)
            best_score = None
            best_fields = None
            for line in fin:
                data = json.loads(line)
                score = data['score']
                fields = data['bump_fields']

                if best_score == None or score < best_score:
                    best_score = score
                    best_fields = fields
            field_list.append(best_fields)

        return field_list


    def hardcoded_bump_fields(self):
        overrides = [{
            '__h_bump_1_field__': 0.0, 
            '__v_bump_1_field__': -0.0, 
            '__h_bump_2_field__': 0.0, 
            '__v_bump_2_field__': -0.0, 
            '__h_bump_3_field__': 0.0, 
            '__v_bump_3_field__': 0.0, 
            '__h_bump_4_field__': -0.0, 
            '__v_bump_4_field__': -0.0, 
            '__h_bump_5_field__': -0.0, 
            '__v_bump_5_field__': 0.0,

        }, {
            "__h_bump_1_field__" : -2.111626439793568e-05 ,
            "__h_bump_2_field__" : 0.05890459353709976 ,
            "__h_bump_3_field__" : 0.0 ,
            "__h_bump_4_field__" : 0.0 ,
            "__h_bump_5_field__" : 0.0 ,
            "__v_bump_1_field__" : -4.564320210898032e-05 ,
            "__v_bump_2_field__" : 1.6563535341074243e-05 ,
            "__v_bump_3_field__" : 0.0 ,
            "__v_bump_4_field__" : 0.0 ,
            "__v_bump_5_field__" : 0.0 ,
        }
        , ]
        return overrides



if __name__ == "__main__":
    Config()
