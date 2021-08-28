import math
import copy
import os
import numpy
import xboa.common
from . import config_double_triplet_baseline as config

class Grid(object):
    def __init__(self, n_points_list, start_list, step_list):
        self.n_points_list = n_points_list
        self.start_list = start_list
        self.step_list = step_list
        self.range_list = range(len(n_points_list))
        self.grid = []
        self.build_grid()

    def greater_than(self, point1, point2):
        for i in self.range_list:
            if point1[i] > point2[i]:
                return True
            elif point1[i] < point2[i]:
                return False
        # equal
        return False

    def add_one(self, point):
        for i in reversed(self.range_list):
            if point[i]+1 < self.n_points_list[i]:
                point[i] += 1
                return point
            else:
                point[i] = 0
        point[0] = self.n_points_list[0]
        return point

    def build_grid(self):
        self.grid = []
        point = [0 for i in self.range_list]
        point[-1] = -1
        max_point = [i-1 for i in self.n_points_list]
        while self.greater_than(max_point, point):
            self.add_one(point)
            self.grid.append(copy.deepcopy(point))
        for item in self.grid:
            for i in self.range_list:
                item[i] = self.start_list[i]+self.step_list[i]*item[i]

    def __getitem__(self, index):
        return self.grid[index]    

def get_seed_fields_grid(n_points_list, start_list, step_list):
    defaults = {
        "__h_bump_1_field__" : 0.0 ,
        "__h_bump_2_field__" : 0.0 ,
        "__h_bump_3_field__" : 0.0 ,
        "__h_bump_4_field__" : 0.0 ,
        "__h_bump_5_field__" : 0.0 ,
        "__v_bump_1_field__" : 0.0 ,
        "__v_bump_2_field__" : 0.0 ,
        "__v_bump_3_field__" : 0.0 ,
        "__v_bump_4_field__" : 0.0 ,
        "__v_bump_5_field__" : 0.0 ,
    }
    targets = ["__h_bump_1_field__", "__h_bump_2_field__", "__v_bump_1_field__", "__v_bump_2_field__"]
    grid = Grid(n_points_list, start_list, step_list)
    seed_fields = []
    for item in grid:
        new_dict = dict(zip(targets, item))
        seed_fields.append(copy.deepcopy(defaults))
        seed_fields[-1].update(new_dict)
    return seed_fields

def get_seed_fields_rotation():
    defaults = {
        "__h_bump_1_field__" : 0.0 ,
        "__h_bump_2_field__" : 0.0 ,
        "__h_bump_3_field__" : 0.0 ,
        "__h_bump_4_field__" : 0.0 ,
        "__h_bump_5_field__" : 0.0 ,
        "__v_bump_1_field__" : 0.0 ,
        "__v_bump_2_field__" : 0.0 ,
        "__v_bump_3_field__" : 0.0 ,
        "__v_bump_4_field__" : 0.0 ,
        "__v_bump_5_field__" : 0.0 ,
    }
    target_1 = "__h_bump_1_field__"
    target_2 = "__v_bump_1_field__"
    seed_fields = []
    for i in range(361):
        fields = copy.deepcopy(defaults)
        phi = math.radians(i)
        fields[target_1] = 0.001*math.sin(phi)
        fields[target_2] = 0.001*math.cos(phi)
        seed_fields.append(fields)
    return seed_fields


class Config(config.Config):
    def __init__(self):
        super(Config, self).__init__()
        co = [3778.0502890298144, 0.00026166947527794804, -103.82702702667791, -5.677699781658774e-05]
        
        ring_tof = 1016.091866
        n_turns = 5
        seed_fields = get_seed_fields_grid([1, 1, 1, 6], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.01])+\
                      get_seed_fields_grid([1, 1, 1, 5], [0.0, 0.0, 0.0, -0.01], [0.0, 0.0, 0.0, -0.01])
        for row in seed_fields:
            for key in sorted(row.keys(), key=lambda key: [x for x in reversed(key)]):
                print(key+format(row[key], "5.3g"), end="")
            print()
        self.run_control["output_dir"] = os.path.join(os.getcwd(), "output/double_triplet_baseline/single_turn_injection/find_bump_parameters_12")
        self.run_control["find_closed_orbits_4d"] = True
        self.run_control["find_da"] = False
        self.run_control["find_bump_parameters"] = True
        self.run_control["track_bump"] = False
        self.tracking["ignore_events"] = []
        self.find_closed_orbits["seed"] = [co]
        self.find_closed_orbits["output_file"] = "closed_orbits_cache"
        self.find_closed_orbits["probe_files"] = ["RINGPROBE*.h5"]
        self.find_closed_orbits["max_iterations"] = 4
        self.find_closed_orbits["us_cell"] = 0
        self.find_closed_orbits["ds_cell"] = 1
        self.find_closed_orbits["subs_overrides"]["__n_turns__"] = 0.5
        self.find_closed_orbits["subs_overrides"]["__do_bump__"] = True
        self.find_closed_orbits["final_subs_overrides"]["__do_bump__"] = True
        self.find_closed_orbits["final_subs_overrides"]["__n_turns__"] = 1
        self.find_closed_orbits["final_subs_overrides"]["__do_magnet_field_maps__"] = True
        self.find_closed_orbits["final_subs_overrides"]["__do_bump__"] = True
        #self.find_closed_orbits["plotting_subs"] = None

        bump = [(0, [0.0*yi, 0.0, 0.0, 0.0]) for yi in range(1)]
        self.find_bump_parameters = {
            "n_h_bumps":5,
            "n_v_bumps":5,
            "amp":[],
            "output_file":"find_bump_parameters",
            "closed_orbit":co,
            "magnet_min_field":-0.2,
            "magnet_max_field":+0.2,
            "max_iterations":1000,
            "field_tolerance":1e-4,
            "amplitude_tolerance":1.,
            "tm_source":"closed_orbits_cache",
            "subs_overrides":{
                "__n_turns__":0.9,
                "__do_magnet_field_maps__":False,
                "__max_x_power__":10,
                "__do_bump__":True,
                "__spt_frequency__":10,
            },
            "final_subs_overrides":None,
            "bump":bump, # beam at foil: 0   3744.07267    3.843    89.83    1.158
            "staged_optimisation":[{ # get foil position
                    "seed_fields":seed,
                    "target_orbit":dict([
                        #(0, [3744.07267, 3.843, -89.83, -1.158]),
                        ]+[(i, co) for i in (1, 2, 7, 8)]),
                    "position_tolerance":0.01,
                    "momentum_tolerance":0.01,
                    "fix_bumps":["__v_bump_1_field__", "__h_bump_1_field__",
                                 "__v_bump_2_field__", "__h_bump_2_field__",
                                 "__v_bump_3_field__", "__h_bump_3_field__",
                                 "__v_bump_-4_field__", "__h_bump_-4_field__",
                                 "__v_bump_-5_field__", "__h_bump_-5_field__",],
                } for seed in seed_fields],
            "target_fields":{},
            "seed_errors":[1e-3]*10,
            "ref_probe_files":["FOILPROBE_1.h5", "RINGPROBE*.h5"], # sorted alphanumerically
            "run_dir":"tmp/find_bump/",
            "energy":3.0,
            "min_time_delta":0., # minimum time between probes
            "target_n_hits":3,
            "penalty_factor":1e9, # penalty = p_f^(number of missed stations)
            "algorithm":"migrad",
            "stop_on_fail":False,
        }
        self.track_bump["input_file"] = "find_bump_parameters_001.out"
        self.track_bump["injection_orbit"] = [3744.0726767615297, 1.6650418416417807, -49.829990504267954, -6.205532738677988]
        self.track_bump["foil_phi"] = 110.8/36.0
        self.track_bump["foil_optimisation_stage"] = 0
        self.track_bump["field_optimisation_stage"] = 0
        self.track_bump["proton_orbit_station"] = 1
        self.track_bump["proton_orbit_phi"]= 0.
        self.track_bump["subs_overrides"]["__magnet_width__"] = 1.0
        self.track_bump["subs_overrides"]["__do_bump__"] = True
        self.track_bump["subs_overrides"]["__septum_field__"] = -4.0
        self.track_bump["subs_overrides"]["__septum_length__"] = 0.1 # m
        self.track_bump["subs_overrides"]["__septum_fringe__"] = 0.1 # fraction of septum_length
        self.track_bump["subs_overrides"]["__septum_width__"] = 0.1 # m
        self.track_bump["subs_overrides"]["__septum_phi__"] = 2.93 # fraction of cell length
        self.track_bump["subs_overrides"]["__septum_dr__"] = 0.0 # m
        self.track_bump["subs_overrides"]["__n_turns__"] = 0.9
        self.track_bump["subs_overrides"]["__septum_dphi__"] = 0.0

    def setup_field_test(self, subs, v_field, h_field, v_phi, h_phi, will_do_main_magnets):      
        overrides = {
            "__n_turns__":0.001,
            "__do_magnet_field_maps__":"False",
            "__do_bump__":"True",
            "__do_rf__":"False",
            "__bump_order__":4,
            "__h_bump_fringe__": 0.1,
            "__v_bump_fringe__": 0.1,
            "__h_bump_length__": 0.2,
            "__v_bump_length__": 0.2,
            "__h_bump_1_field__": h_field,
            "__v_bump_1_field__": v_field, 
            "__h_bump_1_phi__":h_phi,
            "__v_bump_1_phi__":v_phi,
        }
        if not will_do_main_magnets:
            overrides["__bf__"]  = 0.0
            overrides["__bd__"]  = 0.0
        subs.update(overrides)
        return subs

    def setup_field_values(self, subs):
        overrides = {
            '__h_bump_1_field__': 0.051769563021587306, 
            '__v_bump_1_field__': -0.0009659838408628829, 
            '__h_bump_2_field__': 0.016977733266592132, 
            '__v_bump_2_field__': -0.0006105373782377965, 
            '__h_bump_3_field__': 0.0, 
            '__v_bump_3_field__': 0.0, 
            '__h_bump_4_field__': -0.015099536808334202, 
            '__v_bump_4_field__': -0.0206250703423142, 
            '__h_bump_5_field__': -0.035805333570538056, 
            '__v_bump_5_field__': 0.03355348296221017,
        }
        subs.update(overrides)
        return subs

    def setup_linear_ramp(self, subs, n_turns, injection_tof):
        overrides = {
            "__h_bump_1_scales__":[1.0]+[1.0*i/n_turns for i in range(n_turns, -1, -1)]+[0.0],
            "__h_bump_2_scales__":[1.0]+[1.0*i/n_turns for i in range(n_turns, -1, -1)]+[0.0],
            "__h_bump_3_scales__":[1.0]+[1.0*i/n_turns for i in range(n_turns, -1, -1)]+[0.0],
            "__h_bump_4_scales__":[1.0]+[1.0*i/n_turns for i in range(n_turns, -1, -1)]+[0.0],
            "__h_bump_5_scales__":[1.0]+[1.0*i/n_turns for i in range(n_turns, -1, -1)]+[0.0],
            "__v_bump_1_scales__":[1.0]+[1.0*i/n_turns for i in range(n_turns, -1, -1)]+[0.0],
            "__v_bump_2_scales__":[1.0]+[1.0*i/n_turns for i in range(n_turns, -1, -1)]+[0.0],
            "__v_bump_3_scales__":[1.0]+[1.0*i/n_turns for i in range(n_turns, -1, -1)]+[0.0],
            "__v_bump_4_scales__":[1.0]+[1.0*i/n_turns for i in range(n_turns, -1, -1)]+[0.0],
            "__v_bump_5_scales__":[1.0]+[1.0*i/n_turns for i in range(n_turns, -1, -1)]+[0.0],
            "__h_bump_1_times__":[-1e9]+[injection_tof*(i+0.100) for i in range(n_turns+1)]+[1e9],
            "__h_bump_2_times__":[-1e9]+[injection_tof*(i+0.200) for i in range(n_turns+1)]+[1e9],
            "__h_bump_3_times__":[-1e9]+[injection_tof*(i+0.307) for i in range(n_turns+1)]+[1e9],
            "__h_bump_4_times__":[-1e9]+[injection_tof*(i+0.400) for i in range(n_turns+1)]+[1e9],
            "__h_bump_5_times__":[-1e9]+[injection_tof*(i+0.500) for i in range(n_turns+1)]+[1e9],
            "__v_bump_1_times__":[-1e9]+[injection_tof*(i+0.100) for i in range(n_turns+1)]+[1e9],
            "__v_bump_2_times__":[-1e9]+[injection_tof*(i+0.200) for i in range(n_turns+1)]+[1e9],
            "__v_bump_3_times__":[-1e9]+[injection_tof*(i+0.307) for i in range(n_turns+1)]+[1e9],
            "__v_bump_4_times__":[-1e9]+[injection_tof*(i+0.400) for i in range(n_turns+1)]+[1e9],
            "__v_bump_5_times__":[-1e9]+[injection_tof*(i+0.500) for i in range(n_turns+1)]+[1e9],
            "__do_bump__":"True",

        }
        overrides = self.setup_field_values(overrides)
        subs.update(overrides)
        return subs

    def setup_stepped_ramp(self, subs, n_turns, injection_tof):
        step_phase_1 = 0.7499
        step_phase_2 = 0.7501
        steps = [None]*2*n_turns+[0.0, 0.0]
        time = [-1e9]+[None]*2*n_turns+[1e9]
        for i in range(n_turns):
            index = n_turns-i
            steps[2*i] = 1.0*index/n_turns
            steps[2*i+1] = 1.0*index/n_turns
            time[2*i+1] = (i+step_phase_1)*injection_tof
            time[2*i+2] = (i+step_phase_2)*injection_tof
        overrides = {
            "__h_bump_1_scales__":steps,
            "__h_bump_2_scales__":steps,
            "__h_bump_3_scales__":steps,
            "__h_bump_4_scales__":steps,
            "__h_bump_5_scales__":steps,
            "__v_bump_1_scales__":steps,
            "__v_bump_2_scales__":steps,
            "__v_bump_3_scales__":steps,
            "__v_bump_4_scales__":steps,
            "__v_bump_5_scales__":steps,
            "__h_bump_1_times__":time,
            "__h_bump_2_times__":time,
            "__h_bump_3_times__":time,
            "__h_bump_4_times__":time,
            "__h_bump_5_times__":time,
            "__v_bump_1_times__":time,
            "__v_bump_2_times__":time,
            "__v_bump_3_times__":time,
            "__v_bump_4_times__":time,
            "__v_bump_5_times__":time,
            "__do_bump__":"True",

        }
        overrides = self.setup_field_values(overrides)
        subs.update(overrides)
        return subs


    def setup_constant(self, subs, scale_factor):
        steps = [scale_factor, scale_factor]
        time = [-1e9, 1e9]
        overrides = {
            "__h_bump_1_scales__":steps,
            "__h_bump_2_scales__":steps,
            "__h_bump_3_scales__":steps,
            "__h_bump_4_scales__":steps,
            "__h_bump_5_scales__":steps,
            "__v_bump_1_scales__":steps,
            "__v_bump_2_scales__":steps,
            "__v_bump_3_scales__":steps,
            "__v_bump_4_scales__":steps,
            "__v_bump_5_scales__":steps,
            "__h_bump_1_times__":time,
            "__h_bump_2_times__":time,
            "__h_bump_3_times__":time,
            "__h_bump_4_times__":time,
            "__h_bump_5_times__":time,
            "__v_bump_1_times__":time,
            "__v_bump_2_times__":time,
            "__v_bump_3_times__":time,
            "__v_bump_4_times__":time,
            "__v_bump_5_times__":time,
            "__do_bump__":"True",

        }
        overrides = self.setup_field_values(overrides)
        subs.update(overrides)
        return subs

