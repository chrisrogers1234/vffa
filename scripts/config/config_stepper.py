import os
from . import config_double_triplet_baseline as config

class Config(config.Config):
    def __init__(self):
        super(Config, self).__init__()
        self.substitution_list = []
        for i in range(3, 4):
            for stepper in ["RK-4"]:#, "LF-2"]:
                self.substitution_list.append(config.get_baseline_substitution())
                self.substitution_list[-1]["__step_size__"] = 0.1**i
                self.substitution_list[-1]["__stepper__"] = stepper
                self.substitution_list[-1]["__spt_frequency__"] = 1
        self.run_control["output_dir"] = os.path.join(os.getcwd(), "output/double_triplet_baseline/stepper")
