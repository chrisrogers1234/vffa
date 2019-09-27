import os
from . import config_baseline

class Config(config_baseline.Config):
    def __init__(self):
        super(Config, self).__init__()
        self.substitution_list = []
        nsteps=101
        for i in range(nsteps):
            length = 0.125 + 0.1*i/float(nsteps-1)
            self.substitution_list.append(config_baseline.get_baseline_substitution())
            self.substitution_list[-1]["__d_end_length__"] = length
            self.substitution_list[-1]["__f_end_length__"] = length
        print("End length scan", [item["__f_end_length__"] for item in self.substitution_list])
        self.run_control["output_dir"] = os.path.join(os.getcwd(), "output/end_length_scan")
