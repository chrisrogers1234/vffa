import os
import config_baseline

class Config(config_baseline.Config):
    def __init__(self):
        super(Config, self).__init__()
        self.substitution_list = []
        for bb_z in [1.0, 2.0, 4.0]:
            for bb_x in [0.25, 0.5, 1.0, 2.0]:
                self.substitution_list.append(config_baseline.get_baseline_substitution())
                self.substitution_list[-1]["__bb_length__"] = bb_z
                self.substitution_list[-1]["__magnet_width__"] = bb_x
        self.run_control["output_dir"] = os.path.join(os.getcwd(), "output/bb_scan")
