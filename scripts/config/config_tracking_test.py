import math
import os
from . import config_double_triplet_baseline as config

N_TURNS = 10

def time_delay(cell_pos):
    n_turns = N_TURNS
    injection_tof = 1016.091866
    n_cells = 10.0
    td_list = [-1e9, injection_tof*(cell_pos)/n_cells, injection_tof*(n_turns+cell_pos/n_cells), 1e9]
    return td_list

class Config(config.Config):
    def __init__(self):
        super(Config, self).__init__()
        self.run_control["output_dir"] = os.path.join(os.getcwd(), "output/double_triplet_baseline/single_turn_injection/tracking/plotting")
        self.run_control["find_closed_orbits_4d"] = True
        self.run_control["find_da"] = False
        self.run_control["find_bump_parameters"] = False
        self.run_control["track_bump"] = False

        linear_ramp_down = [1.0, 1.0, 0.0, 0.0]
        final_subs = {
            "__n_turns__":0.1, 
            "__do_magnet_field_maps__":"True",
            "__do_bump__":"True", 
            "__do_rf__":"False",
            "__h_bump_1_field__": 0.10239815392945606,
            "__v_bump_1_field__": -0.00023924736234892663, 
            "__h_bump_2_field__": 0.03328208445805614,
            "__v_bump_2_field__": -0.00021702862531158829, 
            "__h_bump_3_field__": 0.0,
            "__v_bump_3_field__": 0.0, 
            "__h_bump_4_field__": -0.07744047791062314,
            "__v_bump_4_field__": -0.016620271881508164, 
            "__h_bump_5_field__": 0.0005000437132594016,
            "__v_bump_5_field__": 0.021246033956753285,
            "__h_bump_1_scales__":linear_ramp_down,
            "__h_bump_2_scales__":linear_ramp_down,
            "__h_bump_3_scales__":linear_ramp_down,
            "__h_bump_4_scales__":linear_ramp_down,
            "__h_bump_5_scales__":linear_ramp_down,
            "__v_bump_1_scales__":linear_ramp_down,
            "__v_bump_2_scales__":linear_ramp_down,
            "__v_bump_3_scales__":linear_ramp_down,
            "__v_bump_4_scales__":linear_ramp_down,
            "__v_bump_5_scales__":linear_ramp_down,
            "__h_bump_1_times__":time_delay(1.0),
            "__h_bump_2_times__":time_delay(2.0),
            "__h_bump_3_times__":time_delay(3.07),
            "__h_bump_4_times__":time_delay(4.0),
            "__h_bump_5_times__":time_delay(5.0),
            "__v_bump_1_times__":time_delay(1.0),
            "__v_bump_2_times__":time_delay(2.0),
            "__v_bump_3_times__":time_delay(3.07),
            "__v_bump_4_times__":time_delay(4.0),
            "__v_bump_5_times__":time_delay(5.0),
        }
        self.find_closed_orbits["final_subs_overrides"] = final_subs
        self.find_closed_orbits["plotting_subs"] = final_subs
        self.find_closed_orbits["plotting_subs"]["__do_magnetic_field_maps__"] = False
        self.find_closed_orbits["plotting_subs"]["__hdf5__"] = False

