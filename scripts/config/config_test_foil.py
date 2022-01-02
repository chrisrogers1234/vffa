import copy
import json
import glob
import math
import os
from . import config_arctan_triplet_baseline as config

class Config(config.Config):
    def __init__(self, angle = None, r0 = None):
        """
        Tracking beam to try to understand the cause of reduced DA. Aim is to
        plot tune and Maxwellianness vs u amplitude
        """
        super(Config, self).__init__()
        if angle == None:
            angle = 0
        if r0 == None:
            r0 = 0
        src_dir = os.path.join(os.getcwd(), "output/arctan_baseline/test_foil/")
        self.run_control["output_dir"] = src_dir.replace("find", "track_")
        ring_tof = 1016.091866
        self.substitution_list = [config.get_baseline_substitution(0.22, ring_tof)] # used for plotting unperturbed CO
        self.run_control["find_closed_orbits_4d"] = False
        self.run_control["find_da"] = False
        self.run_control["find_bump_parameters"] = False
        self.run_control["track_bump"] = False
        self.run_control["track_beam"] = True
        self.find_closed_orbits["subs_overrides"]["__n_turns__"] = 0.11
        self.find_closed_orbits["max_iterations"] = 0
        self.find_closed_orbits["do_minuit"] = True
        self.find_closed_orbits["minuit_iterations"] = 10      
        self.find_closed_orbits["us_cell"] = 0
        self.find_closed_orbits["ds_cell"] = 1

        max_amp = 0.01
        self.track_beam = {
            "run_dir":"tmp/track_beam/",
            "save_dir":"track_beam",
            "print_events":[0, 1, -1],
            "settings":[{
                "name":"forwards",
                "direction":"forwards",
                "probe_files":"RINGPROBE01.h5",          
                "beam":{
                    "type":"beam_gen",
                    "closed_orbit_file":"closed_orbits_cache",
                    "eigen_emittances":[[0, 0]]*2+[[max_amp, max_amp]],
                    "n_per_dimension":4,
                    "variables":["x","x'","y","y'"],
                    "amplitude_dist":"uniform", #"grid", # 
                    "phase_dist":"uniform", #"grid", # 
                    "max_amplitude_4d":max_amp, # amplitude_dist != grid
                    "energy":3.0,
                },
                "subs_overrides":{
                    "__n_turns__":0.1,
                    "__hdf5__":"True",
                    "__do_magnet_field_maps__":"False",
                    "__do_bump__":"True",
                    "__step_size__":0.01,
                    "__do_foil__":True,
                },
            },],
        }


        self.track_beam = {
            "run_dir":"tmp/track_beam/",
            "save_dir":"foil_test",
            "print_events":[i for i in range(1)],
            "variables":["x", "x'", "y", "y'", "t", "energy"],
            "settings":[{
                "name":"grid",
                "direction":"forwards",
                "probe_files":"RINGPROBE01.h5",          
                "beam":{
                    "type":"grid",
                    "energy":3.0,
                    "start":[4357.646683446333, 0.0, -116.7090485272821, 0.0, 0.0, 941.272], #[4370.0, 0.0, -125.0, 0.0, 0.0, 941.272],
                    "stop":[4357.646683446333, 0.0, -116.7090485272821, 0.0, 0.0, 941.272],
                    "nsteps":[1, 1, 1, 1, 1001, 1],
                },
                "subs_overrides":{
                    "__n_turns__":0.51,
                    "__hdf5__":"True",
                    "__do_magnet_field_maps__":"False",
                    "__do_bump__":"False",
                    "__do_rf__":True,
                    "__step_size__":0.01,
                    "__do_foil__":True,
                    "__foil_phi__":3.0, #2.905,
                    "__foil_dphi__":0.0, #0.25/2.+0.025*math.pi/10.0,
                    "__foil_min_z__":-0.150,
                    "__foil_max_z__":-0.100,
                    "__foil_inner_radius__":4.30,
                    "__foil_outer_radius__":4.40,
                    "__foil_thickness__":1.0e-7, # 1e-7 metres = 1e-5 cm
                },
            },],
        }


if __name__ == "__main__":
    Config()
