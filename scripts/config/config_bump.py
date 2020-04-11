import math
import os
from . import config_sector_baseline as config

class Config(config.Config):
    def __init__(self):
        super(Config, self).__init__()
        self.substitution_list = []
        self.substitution_list.append(config.get_baseline_substitution())
        self.run_control["output_dir"] = os.path.join(os.getcwd(), "output/sector_baseline/bumps_2")
        self.run_control["find_closed_orbits_4d"] = False
        self.run_control["find_da"] = False
        self.run_control["find_bump_parameters"] = False
        self.run_control["track_bump"] = True

        self.tracking["dt_tolerance"] = 10.
        self.tracking["verbose"] = 0

        tol = 0.01
        co = [3967.6305564128993, -22.599823825512715, 86.41151974330933, -1.3788397649412332]
        bump = []
        for xi in [0, 1, 2, -2, -1]:
            for yi in range(5):
                bump.append((0, [10.0*xi, 0.0, -10.0*yi, 0.0]))
        self.find_bump_parameters["seed_errors"] = [1e-2]*10
        self.find_bump_parameters["closed_orbit"] = co
        self.find_bump_parameters["bump"] = bump
        self.find_bump_parameters["max_iterations"] = 1000
        self.find_bump_parameters["subs_overrides"] = {
            "__n_turns__":0.6,
            "__do_magnet_field_maps__":False,
            "__do_bump__":True,
        }
        self.find_bump_parameters["algorithm"] = "simplex"
        self.find_bump_parameters["staged_optimisation"][0] = { # get foil position
            "seed_fields":{
                        "__h_bump_1_field__":0.0,
                        "__v_bump_1_field__":0.0,
                        "__h_bump_2_field__":0.0,
                        "__v_bump_2_field__":0.0,
                        "__h_bump_3_field__":0.0,
                        "__v_bump_3_field__":0.1,
                        "__h_bump_4_field__":0.0,
                        "__v_bump_4_field__":0.0,
                        "__h_bump_5_field__":0.0,
                        "__v_bump_5_field__":0.0,
            },
            "target_orbit":dict([
                (0, [3898.02318, -17.75, 88.39, 3.935]),
                ]+[(i, co) for i in (1, 2,)]),
            "position_tolerance":tol,
            "momentum_tolerance":tol,
            "fix_bumps":["__v_bump_-1_field__", "__h_bump_-1_field__",
                         "__v_bump_-2_field__", "__h_bump_-2_field__",
                         "__v_bump_3_field__", "__h_bump_3_field__",
                         "__v_bump_4_field__", "__h_bump_4_field__",
                         "__v_bump_5_field__", "__h_bump_5_field__",],
        }
        self.find_bump_parameters["staged_optimisation"][1] = { # recover closed orbit
            "seed_fields":{},
            "target_orbit":dict([ #maps station:target orbit
                ]+[(i, co) for i in (1, 2, 4)]),
            "position_tolerance":tol,
            "momentum_tolerance":tol,
            "fix_bumps":["__v_bump_1_field__", "__h_bump_1_field__",
                         "__v_bump_2_field__", "__h_bump_2_field__",
                         "__v_bump_3_field__", "__h_bump_3_field__",
                         "__v_bump_-4_field__", "__h_bump_-4_field__",
                         "__v_bump_-5_field__", "__h_bump_-5_field__",],
        }

        self.track_bump["input_file"] = "find_bump_parameters_0*.out"
        self.track_bump["injection_orbit"] = [3898.02318, -17.75, 88.39, 3.935]
        self.track_bump["field_optimisation_stage"] = 0
        self.track_bump["field_optimisation_stage"] = 1
        self.track_bump["foil_station"] = 0