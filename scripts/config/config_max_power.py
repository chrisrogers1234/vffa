import os
from . import config_jan_baseline as config

class Config(config.Config):
    def __init__(self):
        super(Config, self).__init__()
        self.substitution_list = []
        for max_power in range(4, 11, 2):
            self.substitution_list.append(config.get_baseline_substitution())
            self.substitution_list[-1]["__step_size__"] = 0.0001
            self.substitution_list[-1]["__max_x_power__"] = 10
            self.substitution_list[-1]["__bb_length__"] = 4.0
            self.substitution_list[-1]["__magnet_width__"] = 2.0
            
            self.substitution_list[-1]["__max_x_power__"] = max_power
        self.run_control["output_dir"] = os.path.join(os.getcwd(), "output/jan/max_x_power")
