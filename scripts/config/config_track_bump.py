import copy
import json
import glob
import math
import os
from . import config_double_triplet_baseline as config

class Config(config.Config):
    def __init__(self):
        ring_tof = 1016.091866
        super(Config, self).__init__()
        self.run_control["output_dir"] = os.path.join(os.getcwd(), "output/double_triplet_baseline/single_turn_injection/bump_14/")
        self.run_control["find_closed_orbits_4d"] = False
        self.run_control["find_da"] = False
        self.run_control["find_bump_parameters"] = False
        self.run_control["track_bump"] = True
        self.find_closed_orbits["final_subs_overrides"]["__n_turns__"] = 10.1
        self.find_closed_orbits["subs_overrides"]["__n_turns__"] = 1.1
        self.find_closed_orbits["max_iterations"] = 5
        self.find_closed_orbits["us_cell"] = 0
        self.find_closed_orbits["ds_cell"] = 8
        self.find_closed_orbits["seed"] = [[3778.1261886709094, 1.1313917415378683e-05, -103.52060453033312, 9.491741235656548e-06]]
        self.find_da["probe_files"] = "RINGPROBE*.h5"
        self.find_da["x_seed"] = 0.1
        self.find_da["y_seed"] = 0.1
        self.find_da["min_delta"] = 0.0001
        self.find_da["required_n_hits"] = 100
        self.tracking["verbose"] = 0


        self.substitution_list = [config.get_baseline_substitution()]
        bump_fields = self.get_bump_subs(-2e6, -1e6)
        self.substitution_list[-1].update(bump_fields)
        self.find_closed_orbits["ds_cell"] = 6

    def get_bump_subs(self, start_time, end_time):
        subs = {"__do_bump__":True}
        bump_fields = self.get_bump_fields()
        time_step = (end_time-start_time)/len(bump_fields)
        for bump_field in reversed(bump_fields):
            # list of field values __h_bump_1_field__
            # list of scalings     __h_bump_1_scales__ 
            # list of times        __h_bump_1_times__
            for key in bump_field.keys():
                scale_key = key.replace("field", "scales")
                if scale_key not in subs:
                    subs[scale_key] = []
                subs[scale_key].append(bump_field[key])
        times = [-1e9]+[start_time+time_step*i for i in range(len(bump_fields)+1)]+[1e9]
        for key in bump_field[key]:
            time_key = key.replace("field", "times")
            subs[time_key] = times
            subs[key] = 1.0
        print(json.dumps(subs))
        return subs


    def get_bump_fields(self):
        file_name_glob = "output/double_triplet_baseline/single_turn_injection/find_bump_parameters_12/find_bump_parameters*.out"
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

if __name__ == "__main__":
    Config()
