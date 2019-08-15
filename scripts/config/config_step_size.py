import os
import config_baseline

class Config(config_baseline.Config):
    def __init__(self):
        super(Config, self).__init__()
        self.substitution_list = []
        for step in [0.01, 0.005, 0.001, 0.0005, 0.0001]:
            self.substitution_list.append(config_baseline.get_baseline_substitution())
            self.substitution_list[-1]["__step_size__"] = step
        self.run_control["output_dir"] = os.path.join(os.getcwd(), "output/step_size_scan_3")
        self.find_closed_orbits["deltas"] = [0.1, 0.001, 0.1, 0.001]

