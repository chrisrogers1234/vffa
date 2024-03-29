import math
import copy
import os
from . import config_double_triplet_baseline as config

class Config(config.Config):
    def __init__(self):
        super(Config, self).__init__()
        sub = config.get_baseline_substitution()
        sub["__h_bump_1_phi__"] = 1.00
        sub["__h_bump_2_phi__"] = 2.00
        sub["__h_bump_3_phi__"] = 3.07
        sub["__h_bump_4_phi__"] = 4.00
        sub["__h_bump_5_phi__"] = 5.00
        sub["__v_bump_1_phi__"] = 1.00
        sub["__v_bump_2_phi__"] = 2.00
        sub["__v_bump_3_phi__"] = 3.07
        sub["__v_bump_4_phi__"] = 4.00
        sub["__v_bump_5_phi__"] = 5.00
        sub["__h_bump_1_dphi__"] = 0.0
        sub["__h_bump_2_dphi__"] = 0.0
        sub["__h_bump_3_dphi__"] = 0.25
        sub["__h_bump_4_dphi__"] = 0.0
        sub["__h_bump_5_dphi__"] = 0.0
        sub["__v_bump_1_dphi__"] = 0.0
        sub["__v_bump_2_dphi__"] = 0.0
        sub["__v_bump_3_dphi__"] = 0.25
        sub["__v_bump_4_dphi__"] = 0.0
        sub["__v_bump_5_dphi__"] = 0.0
        sub["__foil_probe_phi__"] = 3.095
        sub["__foil_probe_dphi__"] = 0.25/2.+0.025*math.pi/10.
        self.substitution_list = [sub]

        self.run_control["output_dir"] = os.path.join(os.getcwd(), "output/double_triplet_baseline/correlated_painting/fields__0.1_T/")
        self.run_control["find_closed_orbits_4d"] = False
        self.run_control["find_da"] = False
        self.run_control["find_bump_parameters"] = True
        self.run_control["track_bump"] = True

        self.tracking["dt_tolerance"] = 1.
        self.tracking["verbose"] = 0
        self.find_closed_orbits["final_subs_overrides"]["__do_bump__"] = True
        self.find_closed_orbits["final_subs_overrides"]["__cartesian_x_min__"] = -5.0
        self.find_closed_orbits["final_subs_overrides"]["__cartesian_dx__"] = 0.01
        self.find_closed_orbits["final_subs_overrides"]["__cartesian_y_min__"] = -5.0
        self.find_closed_orbits["final_subs_overrides"]["__cartesian_dy__"] = 0.01
        self.find_closed_orbits["final_subs_overrides"]["__n_turns__"] = 1.1

        co = self.find_closed_orbits["seed"][0]
        bumps = BUMP_GLOBAL
        momentum = 75.0
        delta = [-0.0, 0.52, 0.0, -1.13]
        bumps = [bumps[0],]# bumps[12], bumps[24], bumps[32]]
        bumps = [[a_bump[0], a_bump[1]*momentum, a_bump[2], a_bump[3]*momentum] for a_bump in bumps]
        bumps = [(0, [a_bump[i]-bumps[-1][i]+delta[i] for i in range(4)]) for a_bump in bumps] 
        self.find_bump_parameters = {
            "n_h_bumps":5,
            "n_v_bumps":5,
            "amp":[],
            "output_file":"find_bump_parameters",
            "closed_orbit":co,
            "magnet_min_field":-1.0,
            "magnet_max_field":+1.0,
            "max_iterations":5000,
            "field_tolerance":1e-4,
            "amplitude_tolerance":1.,
            "tm_source":"closed_orbits_cache",
            "subs_overrides":{
                "__n_turns__":0.9,
                "__do_magnet_field_maps__":False,
                "__max_x_power__":10,
                "__do_bump__":True,
            },
            "final_subs_overrides":{
                "__n_turns__":0.9,
                "__max_x_power__":10,
                "__do_magnet_field_maps__":True,
                "__do_bump__":True,
            },
            "bump":bumps, # beam at foil: 0   3744.07267    3.843    89.83    1.158
            "staged_optimisation":[{ # get foil position
                    "seed_fields":{
                        "__h_bump_1_field__": 0.05882128321966218, 
                        "__v_bump_1_field__": -0.029353404772833436, 
                        "__h_bump_2_field__": 0.12655412972527147, 
                        "__v_bump_2_field__": 0.03284894560356055, 
                        "__h_bump_3_field__": 0.0, 
                        "__v_bump_3_field__": -0.15, 
                        "__h_bump_4_field__": 0.1456204266846064, 
                        "__v_bump_4_field__": -0.003044975951230655, 
                        "__h_bump_5_field__": -0.04394816758984266, 
                        "__v_bump_5_field__": 0.015611821519977864
                    },
                    "target_orbit":dict([
                        (0, [3744.07267, 3.843, -89.83, -1.158]),
                        ]+[(i, co) for i in (1,)]),
                    "position_tolerance":0.01,
                    "momentum_tolerance":0.01,
                    "fix_bumps":["__v_bump_-1_field__", "__h_bump_-1_field__",
                                 "__v_bump_-2_field__", "__h_bump_-2_field__",
                                 "__v_bump_3_field__", "__h_bump_3_field__",
                                 "__v_bump_4_field__", "__h_bump_4_field__",
                                 "__v_bump_5_field__", "__h_bump_5_field__",],
                }, { # for this optimisation, just want to move the beam using kicker and catch the beam using kickers
                    "seed_fields":{},
                    "target_orbit":dict([ #maps station:target orbit
                        ]+[(i, co) for i in (1, 7, 8)]),
                    "position_tolerance":0.01,
                    "momentum_tolerance":0.001,
                    "fix_bumps":["__v_bump_1_field__", "__h_bump_1_field__",
                                 "__v_bump_2_field__", "__h_bump_2_field__",
                                 "__v_bump_3_field__", "__h_bump_3_field__",
                                 "__v_bump_-4_field__", "__h_bump_-4_field__",
                                 "__v_bump_-5_field__", "__h_bump_-5_field__",],
                },
            ],
            "target_fields":{},
            "seed_errors":[1e-3]*10,
            "ref_probe_files":["FOILPROBE.h5", "RINGPROBE*.h5"], # sorted alphanumerically
            "run_dir":"tmp/find_bump/",
            "energy":3.0,
            "min_time_delta":0., # minimum time between probes
            "target_n_hits":2,
            "penalty_factor":1e9, # penalty = p_f^(number of missed stations)
            "algorithm":"migrad",
        }
        self.track_bump["input_file"] = "find_bump_parameters_00?.out"
        self.track_bump["injection_orbit"] = [3727.076958276747, 2.784588081668254, -44.251054473980474, -2.0746411415980455]
        self.track_bump["foil_phi"] = 110.8/36.0
        self.track_bump["foil_optimisation_stage"] = 0
        self.track_bump["field_optimisation_stage"] = 1
        self.track_bump["proton_orbit_station"] = 1
        self.track_bump["proton_orbit_phi"]= 0.
        self.track_bump["subs_overrides"]["__magnet_width__"] = 1.0
        self.track_bump["subs_overrides"]["__do_bump__"] = True

BUMP_GLOBAL = [
    [
        0.0,
        0.0,
        0.0,
        0.0
    ],
    [
        3.469234850148379,
        0.004295914417639943,
        -6.690975935641678,
        -0.0005806951599070954
    ],
    [
        4.90623897613723,
        0.006075340432220524,
        -9.462468913696469,
        -0.0008212269707450272
    ],
    [
        6.008891023845594,
        0.007440742036320049,
        -11.589110272752094,
        -0.001005793520668423
    ],
    [
        6.938469700296758,
        0.008591828835279886,
        -13.381951871283356,
        -0.0011613903198141907
    ],
    [
        7.7574449548430735,
        0.009605956663364335,
        -14.961477027910048,
        -0.0012984738517573756
    ],
    [
        8.497855180744393,
        0.010522798301883412,
        -16.38947692356337,
        -0.0014224068378762677
    ],
    [
        9.178732653171041,
        0.01136592120269216,
        -17.702658354025594,
        -0.0015363749806530597
    ],
    [
        9.81247795227446,
        0.012150680864441049,
        -18.924937827392938,
        -0.0016424539414900544
    ],
    [
        10.407704550445137,
        0.012887743252919831,
        -20.072927806925033,
        -0.001742085479721286
    ],
    [
        10.970683864501813,
        0.013584874192898045,
        -21.158723726003895,
        -0.001836319331542112
    ],
    [
        11.506150307567015,
        0.014247936254789815,
        -22.191456659820602,
        -0.001925947962987265
    ],
    [
        12.017782047691188,
        0.014881484072640097,
        -23.17822054550419,
        -0.002011587041336846
    ],
    [
        12.50850413883661,
        0.015489139707805837,
        -24.124656818851708,
        -0.0020937261744587926
    ],
    [
        12.980688203511273,
        0.016073839913711172,
        -25.035339534320368,
        -0.0021727623345302585
    ],
    [
        13.436288798707055,
        0.01663800499625183,
        -25.914038368614804,
        -0.0022490226835434326
    ],
    [
        13.876939400593516,
        0.017183657670559772,
        -26.76390374256671,
        -0.0023227806396283815
    ],
    [
        14.304021727235625,
        0.017712508902543268,
        -27.58760052111659,
        -0.0023942674805818914
    ],
    [
        14.718716928411691,
        0.018226021296661572,
        -28.38740674108941,
        -0.002463680912235083
    ],
    [
        15.122044123206257,
        0.0187254568165919,
        -29.165287937124578,
        -0.0025311915190382207
    ],
    [
        15.514889909686147,
        0.01921191332672867,
        -29.922954055820096,
        -0.0025969477035147513
    ],
    [
        15.898031304383732,
        0.01968635299788718,
        -30.661903698205965,
        -0.0026610795259687505
    ],
    [
        16.27215381566463,
        0.02014962468735108,
        -31.383458977133042,
        -0.0027237017296814266
    ],
    [
        16.63786585611668,
        0.020602481785571242,
        -32.088793313877176,
        -0.002784916153317569
    ],
    [
        16.995710361488786,
        0.021045596603766824,
        -32.77895384712674,
        -0.0028448136757525354
    ],
    [
        16.995710361488786,
        0.021045596603766824,
        -34.378953847126745,
        -0.0028448136757525354
    ],
    [
        16.995710361488786,
        0.021045596603766824,
        -35.978953847126746,
        -0.0028448136757525354
    ],
    [
        16.995710361488786,
        0.021045596603766824,
        -37.57895384712675,
        -0.0028448136757525354
    ],
    [
        16.995710361488786,
        0.021045596603766824,
        -39.17895384712675,
        -0.0028448136757525354
    ],
    [
        16.995710361488786,
        0.021045596603766824,
        -40.77895384712675,
        -0.0028448136757525354
    ],
    [
        16.995710361488786,
        0.021045596603766824,
        -42.37895384712675,
        -0.0028448136757525354
    ],
    [
        16.995710361488786,
        0.021045596603766824,
        -43.97895384712675,
        -0.0028448136757525354
    ],
    [
        16.995710361488786,
        0.021045596603766824,
        -45.578953847126755,
        -0.0028448136757525354
    ],
    [
        16.995710361488786,
        0.021045596603766824,
        -47.178953847126756,
        -0.0028448136757525354
    ],
    [
        16.995710361488786,
        0.021045596603766824,
        -48.77895384712676,
        -0.0028448136757525354
    ],
    [
        16.995710361488786,
        0.021045596603766824,
        -50.37895384712676,
        -0.0028448136757525354
    ],
    [
        16.995710361488786,
        0.021045596603766824,
        -51.97895384712676,
        -0.0028448136757525354
    ],
    [
        16.995710361488786,
        0.021045596603766824,
        -53.57895384712676,
        -0.0028448136757525354
    ],
    [
        16.995710361488786,
        0.021045596603766824,
        -55.17895384712676,
        -0.0028448136757525354
    ],
    [
        16.995710361488786,
        0.021045596603766824,
        -56.778953847126765,
        -0.0028448136757525354
    ],
    [
        16.995710361488786,
        0.021045596603766824,
        -58.378953847126766,
        -0.0028448136757525354
    ],
    [
        16.995710361488786,
        0.021045596603766824,
        -59.97895384712677,
        -0.0028448136757525354
    ],
    [
        16.995710361488786,
        0.021045596603766824,
        -61.57895384712677,
        -0.0028448136757525354
    ],
    [
        16.995710361488786,
        0.021045596603766824,
        -63.17895384712677,
        -0.0028448136757525354
    ],
    [
        16.995710361488786,
        0.021045596603766824,
        -64.77895384712677,
        -0.0028448136757525354
    ],
    [
        16.995710361488786,
        0.021045596603766824,
        -66.37895384712677,
        -0.0028448136757525354
    ],
    [
        16.995710361488786,
        0.021045596603766824,
        -67.97895384712676,
        -0.0028448136757525354
    ],
    [
        16.995710361488786,
        0.021045596603766824,
        -69.57895384712675,
        -0.0028448136757525354
    ],
    [
        16.995710361488786,
        0.021045596603766824,
        -71.17895384712675,
        -0.0028448136757525354
    ],
    [
        16.995710361488786,
        0.021045596603766824,
        -72.77895384712674,
        -0.0028448136757525354
    ]
]
