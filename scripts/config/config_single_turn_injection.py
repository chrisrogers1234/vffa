import copy
import json
import glob
import math
import os
from . import config_arctan_triplet_baseline as config

class Config(config.Config):
    def __init__(self, max_r0):
        """
        Tracking beam to try to understand the cause of reduced DA. Aim is to
        plot tune and Maxwellianness vs u amplitude
        """
        super(Config, self).__init__()
        src_dir = os.path.join(os.getcwd(), "output/arctan_baseline/single_turn_injection/tracking_simulation_full/"+max_r0+"-mm")
        self.run_control["output_dir"] = src_dir.replace("find", "track_")
        self.substitution_list = [config.get_baseline_substitution()] # used for plotting unperturbed CO
        injection_tof = self.substitution_list[-1]["__injection_tof__"]

        file_list = [
            "output/arctan_baseline/single_turn_injection/bump_quest_v10/find_bump_r_"+str(ri*-5)+"_theta_90/find_bump_parameters_001.out" \
            for ri in range(13) if abs(ri*-5) <= float(max_r0)
        ]
        file_list = [f for f in reversed(self.sort_list(file_list))]
        for a_file in file_list:
            print(a_file)
        turn_list = [i*2.5 for i, a_file in enumerate(file_list)]
        self.ramp_fields(
                self.substitution_list[0],
                file_list,
                turn_list,
                self.substitution_list[-1]["__injection_tof__"],
                will_step = False,
            )
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

        co = self.find_closed_orbits["seed"][0]
        p1 = (941.272**2-938.272**2)**0.5*(1.0+0.0013)
        p0 = (941.272**2-938.272**2)**0.5*(1.0-0.0013)
        e1 = (p1**2+938.272**2)**0.5
        e0 = (p0**2+938.272**2)**0.5
        n_turns = max(turn_list)+100
        print("Tracking e0", e0, "to e1", e1)
        self.track_beam = {
            "run_dir":"tmp/track_beam/",
            "save_dir":"rf-and-foil",
            "print_events":"all", #[i for i in range(4)],
            "variables":["x", "x'", "y", "y'", "t", "energy"],
            "settings":[{
                "name":"grid",
                "direction":"forwards",
                "probe_files":"RINGPROBE04.h5",        
                "beam":{
                    "type":"grid",
                    "energy":3.0,
                    # phase 0.5 to get in the middle of the (second) RF bucket - 0.3 because we start at cell 3/10
                    "start":[co[0], 0.0, co[2]+0.0, 0.0, injection_tof*0.5-150.0, e0], #[4370.0, 0.0, -125.0, 0.0, 0.0, 941.272],
                    "stop":[co[0], 0.0, co[2]+0.0, 0.0, injection_tof*0.5+150.0, e1],
                    "nsteps":[1, 1, 1, 1, 40+1, 40+1],
                    "reference":[co[0], 0.0, co[2]+0.0, 0.0, injection_tof, 941.272], # reference starts one rf period before t=0 so we don't get any unexpected offsets
                },
                "subs_overrides":{
                    "__n_turns__":n_turns,
                    "__hdf5__":True,
                    "__do_magnet_field_maps__":False,
                    "__do_bump__":True,
                    "__do_rf__":True,
                    "__step_size__":0.01,
                    #"__bf__":0.0,
                    #"__bd__":0.0,
                    "__do_foil__":True,
                    "__foil_phi__":3.05, #2.905,
                    "__foil_dphi__":-0.05, #0.25/2.+0.025*math.pi/10.0,
                    "__foil_min_z__":co[2]/1e3+0.0-0.1/2.0,
                    "__foil_max_z__":co[2]/1e3+0.0+0.1/2.0,
                    "__foil_inner_radius__":co[0]/1e3-0.10,
                    "__foil_outer_radius__":co[0]/1e3-0.020,
                    "__foil_thickness__":1.0e-7, # 1e-7 metres = 1e-5 cm
                    "__beam_phi_init__":0.0, # fraction of cell length
                },
            },],
        }


if __name__ == "__main__":
    Config()
