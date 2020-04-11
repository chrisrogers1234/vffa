import os
from . import config_sector_baseline as config

class Config(config.Config):
    def __init__(self):
        super(Config, self).__init__()
        self.substitution_list = []
        bump_fields = {
            "__h_bump_1_field__" : [0.06073687732641453, 0.04379017068809432, 0.028614168269228912, 0.011231241520031077, -0.008199263136377866, -0.02594158133328861] ,
            "__v_bump_1_field__" : [-0.037442781564765326, -0.02833009075163595, -0.017242412445350297, -0.007431392191126207, 0.0006470717160471651, 0.011223414799383491] ,
            "__h_bump_2_field__" : [-0.06472792547390327, -0.05029789477896607, -0.032509335326137645, -0.015225484594031613, 0.0013565205756123078, 0.022979873179991728] ,
            "__v_bump_2_field__" : [-0.045122469555803235, -0.059755630755301326, -0.0744122527655584, -0.08837225801293513, -0.09921804456406491, -0.11128452365064767] ,
            "__h_bump_3_field__" : [0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ,
            "__v_bump_3_field__" : [0.1, 0.1, 0.1, 0.1, 0.1, 0.1] ,
            "__h_bump_4_field__" : [0.09161860338129912, 0.07467568898282195, 0.055299486954297095, 0.03715800974875094, 0.020737984139591603, 0.0008515138828393276] ,
            "__v_bump_4_field__" : [0.3135470160871101, 0.2513269168247323, 0.18874610063143038, 0.12748945233184372, 0.06138416883011821, 0.0004053908773009596] ,
            "__h_bump_5_field__" : [-0.062398789607350125, -0.05083473769685831, -0.038947472860495935, -0.02647635932232717, -0.014442254432458057, -0.0007520071449533816] ,
            "__v_bump_5_field__" : [-0.32210247240223033, -0.2580368464096914, -0.19367391722521354, -0.130213937391785, -0.06304208532560907, -7.617137287319053e-05] ,
        }
        nsteps = len(bump_fields["__h_bump_1_field__"])
        for i in range(4, 5):
            subs = config.get_baseline_substitution()
            for key in bump_fields:
                subs[key] = bump_fields[key][i]
            self.substitution_list.append(subs)
            self.run_control["find_closed_orbits_4d"] = False
            self.run_control["find_tune"] = True
            self.run_control["find_da"] = False
            self.run_control["find_bump_parameters"] = False
            self.find_closed_orbits["deltas"] = [0.1, 0.1, 0.1, 0.1]
            self.find_closed_orbits["subs_overrides"]["__step_size__"] = 0.001
            self.find_closed_orbits["subs_overrides"]["__n_turns__"] = 1.21
            self.find_closed_orbits["subs_overrides"]["__do_bump__"] = True
            self.find_closed_orbits["final_subs_overrides"]["__n_turns__"] = 0.01
            self.find_closed_orbits["final_subs_overrides"]["__do_bump__"] = True
            self.find_closed_orbits["us_cell"] = 1
            self.find_closed_orbits["ds_cell"] = 2
            self.find_closed_orbits["probe_files"] = "RINGPROBE01.h5"
            self.find_closed_orbits["max_iterations"] = 1
            self.find_tune["delta_1"] = 5.
            self.find_tune["delta_2"] = 5.


        self.run_control["output_dir"] = os.path.join(os.getcwd(), "output/sector/track_bump-0.1T")
