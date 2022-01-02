import math
import os
from . import config_arctan_triplet_baseline as config

class Config(config.Config):
    def __init__(self):
        super(Config, self).__init__()
        self.run_control["output_dir"] = os.path.join(os.getcwd(), "output/arctan_baseline/bb_length")
        self.run_control["find_closed_orbits_4d"] = True
        self.run_control["find_da"] = False
        self.run_control["find_bump_parameters"] = False
        self.run_control["track_bump"] = False

        self.substitution_list = []
        self.find_closed_orbits["do_minuit"] = True
        self.find_closed_orbits["max_iterations"] = 0
        self.find_closed_orbits["seed"] = [
            [4355.014753021532, 0.0, -444.0121684246176, 0.0]
        ]
        bb_list = [bbi for bbi in range(3, 5)]
        bb_list = [2*math.sin(math.radians(36))*2.8]+bb_list
        for bb_length in bb_list:
            sub = config.get_baseline_substitution(0.22)
            sub["__bb_length__"] = bb_length
            self.substitution_list.append(sub)
