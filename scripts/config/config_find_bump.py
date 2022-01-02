import math
import copy
import os
from . import config_arctan_triplet_baseline as config

class Config(config.Config):
    def __init__(self, angle = None, r0 = None):
        super(Config, self).__init__()
        if angle == None:
            angle = 180
        if r0 == None:
            r0 = 20
        ring_tof = 1016.091866
        n_angles = 1
        self.run_control["output_dir"] = os.path.join(os.getcwd(), 
                "output/arctan_baseline/single_turn_injection/bump_quest_v12/find_bump_r_"+str(r0)+"_theta_"+str(angle))
        angle = float(angle)
        r0 = int(r0)
        self.run_control["find_closed_orbits_4d"] = False
        self.run_control["find_da"] = False
        self.run_control["find_bump_parameters"] = True
        self.run_control["track_bump"] = True
        self.find_closed_orbits["max_iterations"] = 2
        self.find_closed_orbits["subs_overrides"]["__n_turns__"] = 0.11
        self.find_closed_orbits["subs_overrides"]["__step_size__"] = 0.01
        self.find_closed_orbits["subs_overrides"]["__do_bump__"] = False
        self.find_closed_orbits["final_subs_overrides"]["__n_turns__"] = 1
        self.find_closed_orbits["final_subs_overrides"]["__step_size__"] = 0.01
        self.find_closed_orbits["final_subs_overrides"]["__do_magnet_field_maps__"] = True
        self.find_closed_orbits["plotting_subs"]["__hdf5__"] = False
        co = self.find_closed_orbits["seed"][0]
        bump = [(4, [i*math.sin(math.radians(angle)), 0.0, i*math.cos(math.radians(angle)), 0.0]) for i in [r0]]
        self.find_bump_parameters = {
            "n_h_bumps":5,
            "n_v_bumps":5,
            "amp":[],
            "output_file":"find_bump_parameters",
            "closed_orbit":co,
            "magnet_min_field":-1.0,
            "magnet_max_field":+1.0,
            "max_iterations":4000,
            "target_score":0.001,
            "field_tolerance":1e-4,
            "amplitude_tolerance":1.,
            "tm_source":"../../..//baseline/closed_orbits_cache",
            "stop_on_fail":False,
            "subs_overrides":{
                "__n_turns__":0.92,
                "__do_magnet_field_maps__":False,
                "__do_bump__":True,
                "__do_foil__":False,
                "__do_rf__":False,
                "__step_size__":0.01,
            },
            "final_subs_overrides":{
                "__n_turns__":0.92,
                "__do_magnet_field_maps__":True,
                "__do_bump__":True,
                "__step_size__":0.01,
            },
            "bump":bump,
            "staged_optimisation":[{ # get foil position
                    "seed_fields": {
                      "__h_bump_1_field__": -6.504114813432604e-05*r0/-25,
                      "__v_bump_1_field__": -0.0017848457494012981*r0/-25,
                      "__h_bump_2_field__": 0.0333635210082841*r0/-25,
                      "__v_bump_2_field__": 0.0031425024028866044*r0/-25,
                      "__h_bump_3_field__": -0.008714704209700774,
                      "__v_bump_3_field__": -0.05011053481841232,
                      "__h_bump_4_field__": 0.03595996696636594*r0/-25,
                      "__v_bump_4_field__": 0.05588885158255952*r0/-25,
                      "__h_bump_5_field__": 0.0427433760001128*r0/-25,
                      "__v_bump_5_field__": -0.046275915940548806*r0/-25
                    },
                    "VERTICAL seed_fields NOT USED": {
                        "__h_bump_1_field__": 0.009389375647487652*r0/-60,
                        "__v_bump_1_field__": 0.020793285608949663*r0/-60,
                        "__h_bump_2_field__": -0.006132646679571363*r0/-60,
                        "__v_bump_2_field__": -0.00872351461761256*r0/-60,
                        "__h_bump_3_field__": -0.0014452685356010075,
                        "__v_bump_3_field__": -0.010060010683301646,
                        "__h_bump_4_field__": -0.00407087907769943*r0/-60,
                        "__v_bump_4_field__": 0.0035875877102686804*r0/-60,
                        "__h_bump_5_field__": 0.008078700510867787*r0/-60,
                        "__v_bump_5_field__": 0.014092378719788634*r0/-60
                    },
                    "zero_fields NOT USED":{
                       '__h_bump_1_field__' : 0.0 ,
                       '__h_bump_2_field__' : 0.0 ,
                       '__h_bump_3_field__' : -0.0 ,
                       '__h_bump_4_field__' : 0.0 ,
                       '__h_bump_5_field__' : -0.0 ,
                       '__v_bump_1_field__' : 0.0 ,
                       '__v_bump_2_field__' : -0.0 ,
                       '__v_bump_3_field__' : -0.0 ,
                       '__v_bump_4_field__' : 0.0 ,
                       '__v_bump_5_field__' : -0.0 ,
                    },
                    "target_orbit":dict([
                        #(0, [3783.04456, 3.96610487, -104.992411, -1.18323871]),
                        ]+[(station, co+tol) for station, tol in (
                            (1, []),
                            #(2, [100, 100, 100, 100]),
                            #(3, [100, 100, 100, 100]),
                            (4, [0.001, 0.001, 0.001, 0.001]), # FOIL
                            #(5, [100, 100, 100, 100]),
                            #(6, [100, 100, 100, 100]),
                            (7, []), # default tolerance
                            (8, []),
                            (9, [])
                        )]),
                    "psv_tolerance":[1e-2, 10, 1e-2, 10], # default if no tolerance is given
                    "fix_bumps":["__v_bump_-1_field__", "__h_bump_-1_field__",
                                 "__v_bump_-2_field__", "__h_bump_-2_field__",
                                 "__v_bump_3_field__", "__h_bump_3_field__",
                                 "__v_bump_-4_field__", "__h_bump_-4_field__",
                                 "__v_bump_-5_field__", "__h_bump_-5_field__",],
                },
            ],
            "target_fields":{},
            "seed_errors":[1e-4]*10,
            "ref_probe_files":["FOILPROBE_1.h5", "RINGPROBE*.h5"], # sorted alphanumerically
            "run_dir":"tmp/find_bump/",
            "energy":3.0,
            "min_time_delta":0., # minimum time between probes
            "target_n_hits":2,
            "penalty_factor":1e9, # penalty = p_f^(number of missed stations)
            "algorithm":"simplex",
        }
        print(bump)
        self.track_bump["input_file"] = "find_bump_parameters_00?.out"
        self.track_bump["injection_orbit"] = None# [3744.0726767615297, 1.6650418416417807, -49.829990504267954, -6.205532738677988]
        self.track_bump["injection_bump"] = 0
        self.track_bump["injection_station"] = 0
        self.track_bump["foil_probe_files"] = ["FOILPROBE*.h5"] # PROBE where the beam is injected
        self.track_bump["foil_phi"] = 110.9/36.0
        self.track_bump["foil_optimisation_stage"] = 0
        self.track_bump["field_optimisation_stage"] = 0
        self.track_bump["proton_orbit_station"] = 1
        self.track_bump["proton_orbit_phi"]= 0.
        self.track_bump["subs_overrides"]["__magnet_width__"] = 1.0
        self.track_bump["subs_overrides"]["__do_bump__"] = True
        self.track_bump["subs_overrides"]["__do_magnet_field_maps__"] = True
        self.track_bump["subs_overrides"]["__septum_field__"] = 0.0
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

