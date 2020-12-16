import math
import copy
import os
from . import config_triplet_baseline as config

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

        # fields: delta 0 0 0 0
        # fields_2: delta -10.0 0 0 0
        # fields_3: delta -10.0 0 0 -1 step 0, 12, 24, 32, 37
        # fields_4: delta -10.0 0 0 -1 step 24, 32, 37
        # trajectory from dp_p_0.0013__col_dens_2e-05__n_tot_50__n_i_25__inj_emit_0.1__k_1.0__tgt_emit_8.0
        self.run_control["output_dir"] = os.path.join(os.getcwd(), "output/triplet_baseline/anticorrelated_painting/fields_3/")
        self.run_control["find_closed_orbits_4d"] = False
        self.run_control["find_da"] = False
        self.run_control["find_bump_parameters"] = True
        self.run_control["track_bump"] = False

        self.tracking["dt_tolerance"] = 1.
        self.tracking["verbose"] = 0
        self.find_closed_orbits["final_subs_overrides"]["__do_bump__"] = True
        self.find_closed_orbits["final_subs_overrides"]["__cartesian_x_min__"] = -5.0
        self.find_closed_orbits["final_subs_overrides"]["__cartesian_dx__"] = 0.01
        self.find_closed_orbits["final_subs_overrides"]["__cartesian_y_min__"] = -5.0
        self.find_closed_orbits["final_subs_overrides"]["__cartesian_dy__"] = 0.01
        self.find_closed_orbits["final_subs_overrides"]["__n_turns__"] = 1.1

        co = self.find_closed_orbits["seed"][0]
        bumps = BUMP_GLOBAL # BUMP_GLOBAL corresponds to setting A_tgt = 8.0 n_i =  25 eps = 0.1
        momentum = 75.0
        delta = [-10.0, 0.66, 0.0, -0.22] #[-0.0, 0.52, 0.0, -1.13]
        bumps = [bumps[0], bumps[12], bumps[25], bumps[32], bumps[37]]
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
            # fields: try with 0.0 delta
            # fields_2: try with delta = [-10.0, 0.0, 0.0, 0.0]
            # <next: try with delta chosen based on v_bump_3 field>
            # beam at foil (0.00 T): 0   3744.07267    3.843    -89.83    -1.158
            # beam at foil (0.10 T):
            # beam at foil (0.15 T):
            # beam at foil (0.20 T):
            "bump":bumps, 
            "staged_optimisation":[{ # get foil position
                    "seed_fields":{
                        "__h_bump_1_field__": -0.12369549075651531, 
                        "__v_bump_1_field__":  -0.05400159877769961, 
                        "__h_bump_2_field__":  0.06507922075334482, 
                        "__v_bump_2_field__":  0.10474410645066845, 
                        "__h_bump_3_field__": 0.0, 
                        "__v_bump_3_field__": -0.1, 
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
            "algorithm":"simplex",
        }
        self.track_bump["input_file"] = "find_bump_parameters_001.out"
        self.track_bump["injection_orbit"] = [3744.0726767615297, 1.6650418416417807, -49.829990504267954, -6.205532738677988]
        self.track_bump["foil_phi"] = 110.8/36.0
        self.track_bump["foil_optimisation_stage"] = 0
        self.track_bump["field_optimisation_stage"] = 0
        self.track_bump["proton_orbit_station"] = 1
        self.track_bump["proton_orbit_phi"]= 0.
        self.track_bump["subs_overrides"]["__magnet_width__"] = 1.0
        self.track_bump["subs_overrides"]["__do_bump__"] = True

BUMP_GLOBAL =  [
        [
            22.67220396322874,
            0.00761296136164872,
            -9.188295597408736,
            -0.001568559847107185
        ],
        [
            21.235716634860687,
            0.0074150848055477515,
            -13.927671370237881,
            -0.0015951368471192632
        ],
        [
            20.35057606105754,
            0.007235701323738961,
            -15.773205778273137,
            -0.0015860731047297315
        ],
        [
            19.54665282424491,
            0.007056172878279011,
            -17.138784786820388,
            -0.0015704888285268319
        ],
        [
            18.778547280181847,
            0.00687447964478903,
            -18.25339842219648,
            -0.0015510987077416412
        ],
        [
            18.028074167783785,
            0.00668963776035996,
            -19.205500846106467,
            -0.0015289127375282517
        ],
        [
            17.285339434765703,
            0.006500952068295453,
            -20.04022673130516,
            -0.0015044096398255827
        ],
        [
            16.54389859259675,
            0.0063078192350183755,
            -20.78415837043503,
            -0.0014778345721823502
        ],
        [
            15.798963806235628,
            0.006109648414110812,
            -21.454376938309103,
            -0.0014493065260858845
        ],
        [
            15.046578097651555,
            0.005905817209624622,
            -22.062491199656574,
            -0.001418864474754858
        ],
        [
            14.283164429532414,
            0.005695638758058923,
            -22.616671365397043,
            -0.0013864886430914611
        ],
        [
            13.505234849089911,
            0.005478330076743567,
            -23.12276610524258,
            -0.0013521096144480244
        ],
        [
            12.70916509860325,
            0.005252975803272236,
            -23.58494865116155,
            -0.0013156101923772873
        ],
        [
            11.890982766945616,
            0.0050184819165424496,
            -24.00609592155957,
            -0.0012768217499859694
        ],
        [
            11.046130172622009,
            0.004773512490412464,
            -24.387998723247023,
            -0.0012355150751530708
        ],
        [
            10.169160287501768,
            0.004516398771613242,
            -24.731447694120448,
            -0.0011913842843452
        ],
        [
            9.253305824927876,
            0.004245002357088674,
            -25.036205901863514,
            -0.0011440204952591223
        ],
        [
            8.289818744444942,
            0.0039564990946820465,
            -25.30084750268544,
            -0.0010928686478756864
        ],
        [
            7.266882163859363,
            0.0036470179182089804,
            -25.522395296411766,
            -0.0010371540844671197
        ],
        [
            6.167672011560703,
            0.003310993141522146,
            -25.695594291849794,
            -0.000975749845635279
        ],
        [
            4.966554201035784,
            0.0029398899506447498,
            -25.811415812743718,
            -0.0009069146460353235
        ],
        [
            3.620574836139668,
            0.002519348832606151,
            -25.853643946361238,
            -0.0008277049633278848
        ],
        [
            2.0462097556555277,
            0.0020213799285294784,
            -25.7894806248088,
            -0.0007323671286257272
        ],
        [
            0.028145658203523105,
            0.0013737343846172202,
            -25.53259973462058,
            -0.000606026938622793
        ],
        [
            -4.698730314345225,
            -0.0001841317378434217,
            -24.165858652936958,
            -0.00029199392974905123
        ],
        [
            -4.698730314345225,
            -0.0001841317378434217,
            -25.76585865293696,
            -0.00029199392974905123
        ],
        [
            -4.698730314345225,
            -0.0001841317378434217,
            -27.36585865293696,
            -0.00029199392974905123
        ],
        [
            -4.698730314345225,
            -0.0001841317378434217,
            -28.965858652936962,
            -0.00029199392974905123
        ],
        [
            -4.698730314345225,
            -0.0001841317378434217,
            -30.565858652936964,
            -0.00029199392974905123
        ],
        [
            -4.698730314345225,
            -0.0001841317378434217,
            -32.16585865293696,
            -0.00029199392974905123
        ],
        [
            -4.698730314345225,
            -0.0001841317378434217,
            -33.76585865293696,
            -0.00029199392974905123
        ],
        [
            -4.698730314345225,
            -0.0001841317378434217,
            -35.365858652936964,
            -0.00029199392974905123
        ],
        [
            -4.698730314345225,
            -0.0001841317378434217,
            -36.965858652936966,
            -0.00029199392974905123
        ],
        [
            -4.698730314345225,
            -0.0001841317378434217,
            -38.56585865293697,
            -0.00029199392974905123
        ],
        [
            -4.698730314345225,
            -0.0001841317378434217,
            -40.16585865293697,
            -0.00029199392974905123
        ],
        [
            -4.698730314345225,
            -0.0001841317378434217,
            -41.76585865293697,
            -0.00029199392974905123
        ],
        [
            -4.698730314345225,
            -0.0001841317378434217,
            -43.36585865293697,
            -0.00029199392974905123
        ],
        [
            -4.698730314345225,
            -0.0001841317378434217,
            -44.96585865293697,
            -0.00029199392974905123
        ],
        [
            -4.698730314345225,
            -0.0001841317378434217,
            -46.565858652936974,
            -0.00029199392974905123
        ],
        [
            -4.698730314345225,
            -0.0001841317378434217,
            -48.165858652936976,
            -0.00029199392974905123
        ],
        [
            -4.698730314345225,
            -0.0001841317378434217,
            -49.76585865293698,
            -0.00029199392974905123
        ],
        [
            -4.698730314345225,
            -0.0001841317378434217,
            -51.36585865293698,
            -0.00029199392974905123
        ],
        [
            -4.698730314345225,
            -0.0001841317378434217,
            -52.96585865293698,
            -0.00029199392974905123
        ],
        [
            -4.698730314345225,
            -0.0001841317378434217,
            -54.56585865293698,
            -0.00029199392974905123
        ],
        [
            -4.698730314345225,
            -0.0001841317378434217,
            -56.16585865293698,
            -0.00029199392974905123
        ],
        [
            -4.698730314345225,
            -0.0001841317378434217,
            -57.765858652936984,
            -0.00029199392974905123
        ],
        [
            -4.698730314345225,
            -0.0001841317378434217,
            -59.365858652936986,
            -0.00029199392974905123
        ],
        [
            -4.698730314345225,
            -0.0001841317378434217,
            -60.96585865293699,
            -0.00029199392974905123
        ],
        [
            -4.698730314345225,
            -0.0001841317378434217,
            -62.56585865293699,
            -0.00029199392974905123
        ],
        [
            -4.698730314345225,
            -0.0001841317378434217,
            -64.16585865293699,
            -0.00029199392974905123
        ]
    ]

