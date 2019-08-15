import os
import config_baseline

class Config(config_baseline.Config):
    def __init__(self):
        super(Config, self).__init__()
        self.substitution_list = []
        for i in range(10):
            self.substitution_list.append(config_baseline.get_baseline_substitution())
            self.substitution_list[-1]["__energy__"] = i*1.0+3.0
        self.run_control["output_dir"] = os.path.join(os.getcwd(), "output/energy_scan")
