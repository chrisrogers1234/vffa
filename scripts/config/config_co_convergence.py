import math
import os
from . import config_double_triplet_baseline as config

class Config(config.Config):
    def __init__(self):
        super(Config, self).__init__()
        delta = 0.001
        self.run_control["output_dir"] = os.path.join(os.getcwd(), "output/double_triplet_baseline/co_convergence_"+str(delta))
        self.run_control["find_closed_orbits_4d"] = True
        self.run_control["find_da"] = False
        self.run_control["find_bump_parameters"] = False
        self.run_control["track_bump"] = False

        self.substitution_list = []
        self.find_closed_orbits["deltas"] = [delta, delta*0.1, delta, delta*0.1]
        self.find_closed_orbits["do_minuit"] = False
        self.find_closed_orbits["minuit_tolerance"] = 1e-14
        self.find_closed_orbits["seed"] = [[3778.0377273690738, 0.0, -103.98094256513187, 0.0]]
        self.find_closed_orbits["max_iterations"] = 1
        for step in [0.1**i for i in range(2, 7)]:
            sub = config.get_baseline_substitution()
            sub["__output_plane_tolerance__"] = 1e-14
            sub["__step_size__"] = step
            self.substitution_list.append(sub)
