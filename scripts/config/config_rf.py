import math
import os
from . import config_double_triplet_baseline as config

class Config(config.Config):   
    def __init__(self):
        super(Config, self).__init__()
        self.run_control["output_dir"] = os.path.join(os.getcwd(), "output/double_triplet_baseline/rf/plotting/")
        self.run_control["find_closed_orbits_4d"] = True
        self.run_control["find_da"] = False
        self.run_control["find_bump_parameters"] = False
        self.run_control["track_bump"] = False
        self.find_closed_orbits["final_subs_overrides"]["__n_turns__"] = 0.001
        self.find_closed_orbits["final_subs_overrides"]["__spt_frequency__"] = 100

        self.substitution_list = []
        print("Running phases")
        sub = config.get_baseline_substitution()
        # 1.0 3.009509181
        # 2.0 3.005872602
        # 3.0 2.999995649
        # 3.5 2.996908856
        # 6.0 2.990505249
        phi0 = 8.0/10.0*math.pi*2.0
        t0 = 1016.091866
        rf_params = self.linear_ramp([0.0, 10e6], [2./t0*1e3, 4./t0*1e3], [phi0, phi0], [1.0e-1, 1.0e-1])
        rf_params = self.load_rf("lattice/2021-02-15_kelliher-rf-programme.txt", 0.1, 2, phi0)
        sub.update(rf_params)
        self.substitution_list.append(sub)

    def linear_ramp(self, time, frequency, phase, voltage):
        rf_params = {
            "__rf_time__":time,
            "__rf_phase__":phase,
            "__rf_efield__":voltage,
            "__rf_frequency__":frequency,
            "__do_rf__":True,
        }
        return rf_params


    def load_rf(self, file_name, rf_length, harmonic_number, phase):
        heading = ["__rf_time__", "__rf_efield__", "__rf_frequency__", "__rf_phase__"]
        fin = open(file_name)
        fin.readline()
        #  turn    V_0 [V] phi_s [deg]    f_rf [Hz]   f_s [Hz] B.A. [Vs] K.E. [eV]

        data = {"__do_rf__":True}
        for head in heading:
            data[head] = []
        time = 0
        for line in fin.readlines():
            words = line.split()
            voltage = float(words[1])
            freq = float(words[3])
            time += (1/freq)*harmonic_number
            data["__rf_time__"].append(time*1e9) # ns
            data["__rf_efield__"].append(voltage*1e-6/rf_length) # MV/m
            data["__rf_phase__"].append(phase) # ns
            data["__rf_frequency__"].append(freq*1e6) # MHz
        return data


