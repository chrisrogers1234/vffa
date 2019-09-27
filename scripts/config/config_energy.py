import os
from . import config_jan_baseline

class Config(config_jan_baseline.Config):
    def __init__(self):
        super(Config, self).__init__()
        self.substitution_list = []
        for i in range(5):
            self.substitution_list.append(config_jan_baseline.get_baseline_substitution())
            self.substitution_list[-1]["__energy__"] = i*100.0+400.0
        self.run_control["output_dir"] = os.path.join(os.getcwd(), "output/jan_energy_scan")
        self.find_closed_orbits["max_iterations"] = 10
        self.run_control["find_closed_orbits_4d"] = False
        self.run_control["find_tune"] = True