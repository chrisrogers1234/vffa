import math
import copy
import os
from . import config_triplet_baseline as config

class Config(config.Config):
    def __init__(self):
        super(Config, self).__init__()
        sub = config.get_baseline_substitution()
        sub["__h_bump_1_phi__"] = 1.00
        sub["__h_bump_2_phi__"] = 2.00
        sub["__h_bump_3_phi__"] = 3.07
        sub["__h_bump_4_phi__"] = 4.00
        sub["__h_bump_5_phi__"] = 5.00
        sub["__v_bump_1_phi__"] = 1.00
        sub["__v_bump_2_phi__"] = 2.00
        sub["__v_bump_3_phi__"] = 3.07
        sub["__v_bump_4_phi__"] = 4.00
        sub["__v_bump_5_phi__"] = 5.00
        sub["__h_bump_1_dphi__"] = 0.0
        sub["__h_bump_2_dphi__"] = 0.0
        sub["__h_bump_3_dphi__"] = 0.25
        sub["__h_bump_4_dphi__"] = 0.0
        sub["__h_bump_5_dphi__"] = 0.0
        sub["__v_bump_1_dphi__"] = 0.0
        sub["__v_bump_2_dphi__"] = 0.0
        sub["__v_bump_3_dphi__"] = 0.25
        sub["__v_bump_4_dphi__"] = 0.0
        sub["__v_bump_5_dphi__"] = 0.0
        sub["__foil_probe_phi__"] = 3.095
        sub["__foil_probe_dphi__"] = 0.25/2.+0.025*math.pi/10.
        self.substitution_list = [sub]

        self.run_control["output_dir"] = os.path.join(os.getcwd(), "output/triplet_baseline/single_turn_injection/septum/")
        self.run_control["find_closed_orbits_4d"] = False
        self.run_control["find_da"] = False
        self.run_control["find_bump_parameters"] = False
        self.run_control["track_bump"] = True

        self.tracking["lattice_file"] = os.path.join(os.getcwd(), "lattice/VerticalTripletFFA1.in")
        self.tracking["dt_tolerance"] = 1.
        self.tracking["verbose"] = 0
        self.find_closed_orbits["final_subs_overrides"]["__do_bump__"] = True
        self.find_closed_orbits["final_subs_overrides"]["__cartesian_x_min__"] = -5.0
        self.find_closed_orbits["final_subs_overrides"]["__cartesian_dx__"] = 0.01
        self.find_closed_orbits["final_subs_overrides"]["__cartesian_y_min__"] = -5.0
        self.find_closed_orbits["final_subs_overrides"]["__cartesian_dy__"] = 0.01
        self.find_closed_orbits["final_subs_overrides"]["__n_turns__"] = 1.1

        co = self.find_closed_orbits["seed"][0]
        bump = [(0, [0.0, 0.0, +40.0*yi, 0.0]) for yi in range(1, 2)]
        self.find_bump_parameters = {
            "n_h_bumps":5,
            "n_v_bumps":5,
            "amp":[],
            "output_file":"find_bump_parameters",
            "closed_orbit":co,
            "magnet_min_field":-1.0,
            "magnet_max_field":+1.0,
            "max_iterations":5000,
            "field_tolerance":1e-4,
            "amplitude_tolerance":1.,
            "tm_source":"closed_orbits_cache",
            "subs_overrides":{
                "__n_turns__":0.9,
                "__do_magnet_field_maps__":False,
                "__max_x_power__":10,
                "__do_bump__":True,
            },
            "final_subs_overrides":{
                "__n_turns__":0.9,
                "__max_x_power__":10,
                "__do_magnet_field_maps__":True,
                "__do_bump__":True,
            },
            "bump":bump, # beam at foil: 0   3744.07267    3.843    89.83    1.158
            "staged_optimisation":[{ # get foil position
                    "seed_fields":{
                        "__h_bump_1_field__": 0.09594564327812449,
                        "__v_bump_1_field__": 0.010822238322893973,
                        "__h_bump_2_field__": 0.04303470161076883,
                        "__v_bump_2_field__": -0.007344463304957749,
                        "__h_bump_3_field__": 0.0,
                        "__v_bump_3_field__": -0.2,
                        "__h_bump_4_field__": 0.08266218278182147,
                        "__v_bump_4_field__": 0.16451116210773065,
                        "__h_bump_5_field__": -0.024656607344269554,
                        "__v_bump_5_field__": -0.09561693485296463
                    },
                    "target_orbit":dict([
                        (0, [3744.07267, 3.843, -89.83, -1.158]),
                        ]+[(i, co) for i in (1, 7, 8, 9)]),
                    "position_tolerance":0.01,
                    "momentum_tolerance":100.0,
                    "fix_bumps":["__v_bump_-1_field__", "__h_bump_-1_field__",
                                 "__v_bump_-2_field__", "__h_bump_-2_field__",
                                 "__v_bump_3_field__", "__h_bump_3_field__",
                                 "__v_bump_-4_field__", "__h_bump_-4_field__",
                                 "__v_bump_-5_field__", "__h_bump_-5_field__",],
                },
            ],
            "void":[{ # for this optimisation, just want to move the beam using kicker and catch the beam using kickers
                    "seed_fields":{},
                    "target_orbit":dict([ #maps station:target orbit
                        ]+[(i, co) for i in (1, 7, 8)]),
                    "position_tolerance":0.01,
                    "momentum_tolerance":0.001,
                    "fix_bumps":["__v_bump_1_field__", "__h_bump_1_field__",
                                 "__v_bump_2_field__", "__h_bump_2_field__",
                                 "__v_bump_3_field__", "__h_bump_3_field__",
                                 "__v_bump_-4_field__", "__h_bump_-4_field__",
                                 "__v_bump_-5_field__", "__h_bump_-5_field__",],
                },
            ],
            "target_fields":{},
            "seed_errors":[1e-3]*10,
            "ref_probe_files":["FOILPROBE.h5", "RINGPROBE*.h5"], # sorted alphanumerically
            "run_dir":"tmp/find_bump/",
            "energy":3.0,
            "min_time_delta":0., # minimum time between probes
            "target_n_hits":4,
            "penalty_factor":1e9, # penalty = p_f^(number of missed stations)
            "algorithm":"migrad",
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


