import math
import os
from . import config_double_triplet_baseline as config

class Config(config.Config):
    def __init__(self):
        super(Config, self).__init__()
        self.run_control["output_dir"] = os.path.join(os.getcwd(), "output/double_triplet_baseline/rf/tracking/plotting")
        self.run_control["find_closed_orbits_4d"] = True
        self.run_control["find_da"] = False
        self.run_control["find_bump_parameters"] = False
        self.run_control["track_bump"] = False
        self.find_closed_orbits["final_subs_overrides"]["__n_turns__"] = 200.1

        self.substitution_list = []
        print("Running phases")
        for phase in [-0.0011294694720536897+i/10.0 for i in range(1)]:
            for voltage in [1.0, 0.1]:
                print(format(phase, "6.2g")+"*2*pi")
                sub = config.get_baseline_substitution()
                sub["__rf_phase__"] = [phase*math.pi*2.0]*4
                sub["__rf_efield__"] = [voltage]*4
                self.substitution_list.append(sub)
#-0.05, 2.96971174
#0, 3.000753706 
