import os
from . import config_jan_baseline as config

class Config(config.Config):
    def __init__(self):
        super(Config, self).__init__()
        self.substitution_list = []
        for step in [0.01, 0.005, 0.001, 0.0005, 0.0001]:
            self.substitution_list.append(config.get_baseline_substitution())
            self.substitution_list[-1]["__step_size__"] = 0.0001
            self.substitution_list[-1]["__max_x_power__"] = 10
            self.substitution_list[-1]["__bb_length__"] = 4.0
            self.substitution_list[-1]["__magnet_width__"] = 2.0
            self.substitution_list[-1]["__step_size__"] = step
        self.run_control["output_dir"] = os.path.join(os.getcwd(), "output/jan/step_size_scan_0.01")
        self.find_closed_orbits["deltas"] = [0.01, 0.01, 0.01, 0.01]

