import math
import os

def get_baseline_substitution(offset):
    field_factor = 1.8
    injection_tof = 1016.091866
    baseline = {
        # ring
        "__n_cells__":10,
        "__radius__":4.04508497187, #m
        "__half_cell_length__":1.25,
        "__n_particles__":3,
        "__drift_radius__":3.74,
        # beam
        "__energy__":3., # MeV
        "__beam_phi_init__":0., # fraction of cell length
        "__beam_charge__":1., # multiples of positron charge
        # tracking
        "__step_size__":0.001, # m
        "__n_turns__":0.2,
        "__hdf5__":True,
        "__stepper__":"RK-4",
        "__spt_frequency__":1,
        # main magnets
        "__bf__":-0.5/field_factor, #-0.25/0.44, # T
        "__bfa__":1.6, #-0.25/0.44, # T
        "__bfb__":-0.6, #-0.25/0.44, # T
        "__bf_offset__":offset,
        "__bd__":0.1/field_factor, #0.25, #+0.20, # T
        "__bda__":1.6, #0.25, #+0.20, # T
        "__bdb__":-0.6, #0.25, #+0.20, # T
        "__bd_offset__":offset,
        "__d_length__":0.24, # m
        "__f_length__":0.40, # m
        "__fd_gap__":0.08, # m # 0.08
        "__drift_length__":1.30, # m # 0.08
        "__f_tilt_angle__":0.0,  # m
        "__m_index__":-1.28, # m^-1 # 1.28
        "__d_end_length__":0.20,  # m
        "__f_end_length__":0.20, # m
        "__max_x_power__":10,
        "__neg_extent__":2.0, # m
        "__pos_extent__":4.0, # m
        # NOTE BB length 0.125 m -> width 1.0 m; length 2.0 m; double for safety
        "__magnet_width__":2.0, # m
        "__bb_length__":8.0, # m
        # field maps
        "__do_magnet_field_maps__":False,
        "__cartesian_x_min__":-5., # m
        "__cartesian_dx__":0.01, # m (1001 steps)
        "__cartesian_y_min__":-5., # m
        "__cartesian_dy__":0.01, # m (1001 steps)
        # rf
        "__do_rf__":False,
        "__rf_phi__":6.0,
        "__rf_length__":0.1,
        "__rf_end_factor__":0.1,
        "__do_rf_fringe_field__":False,
        "__rf_time__":[-1e9, 0.0, 1e6, 1e9],
        "__rf_frequency__":[2.0/injection_tof, 2.0/injection_tof, 4.0/injection_tof, 4.0/injection_tof],
        "__rf_efield__":[1.0]*4,
        "__rf_phase__":[0.0, 0.0, 0.0, 0.0],
        # bump magnets
        "__do_bump__":False,
        "__bump_order__":4,
        "__h_bump_fringe__":0.1,
        "__h_bump_length__":0.2,
        "__v_bump_fringe__":0.1,
        "__v_bump_length__":0.2,
        "__h_bump_1_field__": 0.0,
        "__v_bump_1_field__": 0.0, 
        "__h_bump_2_field__": 0.0,
        "__v_bump_2_field__": 0.0, 
        "__h_bump_3_field__": 0.0, 
        "__v_bump_3_field__": 0.0, 
        "__h_bump_4_field__": 0.0, 
        "__v_bump_4_field__": 0.0, 
        "__h_bump_5_field__": 0.0, 
        "__v_bump_5_field__": 0.0,
        "__h_bump_1_scales__":[1.0, 1.0],
        "__h_bump_2_scales__":[1.0, 1.0],
        "__h_bump_3_scales__":[1.0, 1.0],
        "__h_bump_4_scales__":[1.0, 1.0],
        "__h_bump_5_scales__":[1.0, 1.0],
        "__v_bump_1_scales__":[1.0, 1.0],
        "__v_bump_2_scales__":[1.0, 1.0],
        "__v_bump_3_scales__":[1.0, 1.0],
        "__v_bump_4_scales__":[1.0, 1.0],
        "__v_bump_5_scales__":[1.0, 1.0],
        "__h_bump_1_times__":[-1e9, 1e9],
        "__h_bump_2_times__":[-1e9, 1e9],
        "__h_bump_3_times__":[-1e9, 1e9],
        "__h_bump_4_times__":[-1e9, 1e9],
        "__h_bump_5_times__":[-1e9, 1e9],
        "__v_bump_1_times__":[-1e9, 1e9],
        "__v_bump_2_times__":[-1e9, 1e9],
        "__v_bump_3_times__":[-1e9, 1e9],
        "__v_bump_4_times__":[-1e9, 1e9],
        "__v_bump_5_times__":[-1e9, 1e9],
        "__h_bump_1_phi__":1.00,
        "__h_bump_2_phi__":2.00,
        "__h_bump_3_phi__":3.07,
        "__h_bump_4_phi__":4.00,
        "__h_bump_5_phi__":5.00,
        "__v_bump_1_phi__":1.00,
        "__v_bump_2_phi__":2.00,
        "__v_bump_3_phi__":3.07,
        "__v_bump_4_phi__":4.00,
        "__v_bump_5_phi__":5.00,
        "__h_bump_1_dphi__":0.0,
        "__h_bump_2_dphi__":0.0,
        "__h_bump_3_dphi__":0.25,
        "__h_bump_4_dphi__":0.0,
        "__h_bump_5_dphi__":0.0,
        "__v_bump_1_dphi__":0.0,
        "__v_bump_2_dphi__":0.0,
        "__v_bump_3_dphi__":0.25,
        "__v_bump_4_dphi__":0.0,
        "__v_bump_5_dphi__":0.0,
        "__foil_probe_1_phi__":3.095,
        "__foil_probe_1_dphi__":0.25/2.+0.025*math.pi/10.0,
        "__foil_probe_2_phi__":4.095,
        "__foil_probe_2_dphi__":0.25/2.+0.025*math.pi/10.0,
        "__septum_field__":-0.0,
        "__septum_length__":0.1, # m
        "__septum_fringe__":0.1, # fraction of septum_length
        "__septum_width__":0.5, # m
        "__septum_phi__":0.5, # fraction of cell length
        "__septum_dr__":0.5, # m
    }
    return baseline

class Config(object):
    def __init__(self):
        delta = 0.1
        xp_ratio = 0.1
        tol = 1e-12
        offset = 0.03
        self.find_closed_orbits = {
            "seed":[[3778.1261886709094, 1.1313917415378683e-05, -103.52060453033312, 9.491741235656548e-06]],#, 0., 3.]],
            "deltas":[delta, delta*xp_ratio, delta, delta*xp_ratio],
            "adapt_deltas":False,
            "output_file":"closed_orbits_cache",
            "subs_overrides":{"__n_turns__":0.11, "__do_magnet_field_maps__":"False"},
            "final_subs_overrides":{"__n_turns__":1.001, "__do_magnet_field_maps__":"False",
                    "__do_bump__":"False", "__do_rf__":"False",},
            "plotting_subs":{
                "__n_turns__":0.001, 
                "__do_magnet_field_maps__":"False",
                "__hdf5__":"False"
            },
            "us_cell":0,
            "ds_cell":1,
            "root_batch":0,
            "max_iterations":10, 
            "tolerance":0.00001,
            "do_plot_orbit":False,
            "run_dir":"tmp/find_closed_orbits",
            "probe_files":"RINGPROBE*.h5",
            "overwrite":True,
            "orbit_file":"VerticalSectorFFA-trackOrbit.dat",
            "use_py_tracking":False,
            "py_tracking_phi_list":[0.0, math.pi*2.0/10.0],
        }
        self.find_tune = {
            "run_dir":"tmp/find_tune/",
            "probe_files":"RINGPROBE01.h5",
            "subs_overrides":{"__n_turns__":120.1, "__do_magnet_field_maps__":"True"},
            "root_batch":0,
            "delta_1":5.,
            "delta_2":5.,
            "max_n_cells":0.1,
            "output_file":"find_tune",
            "row_list":None,
            "axis":None,
        }
        self.find_da = {
            "run_dir":"tmp/find_da/",
            "probe_files":"RINGPROBE*.h5",
            "subs_overrides":{"__n_turns__":12.1, "__do_magnet_field_maps__":"False", "__step_size__":0.001},
            "get_output_file":"get_da",
            "scan_output_file":"scan_da",
            "row_list":None,
            "scan_x_list":[],
            "scan_y_list":[],
            "x_seed":1.,
            "y_seed":1.,
            "min_delta":0.9,
            "max_delta":1000.,
            "required_n_hits":100,
            "dt_tolerance":0.5, # fraction of closed orbit dt
            "max_iterations":15,
            "decoupled":True,
        }
        co = self.find_closed_orbits["seed"][0]
        bump = []
        for xi in [0, 1, 2, -2, -1]:
            for yi in range(5):
                bump.append((0, [10.0*xi, 0.0, -10.0*yi, 0.0]))
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
                "__n_turns__":1.1,
                "__do_magnet_field_maps__":False,
                "__max_x_power__":10,
                "__do_bump__":True,
            },
            "final_subs_overrides":{
                "__n_turns__":1.1,
                "__max_x_power__":10,
                "__do_magnet_field_maps__":True,
                "__do_bump__":True,
            },
            "bump":bump,
            "staged_optimisation":[{ # get foil position
                    "seed_fields":{
                                "__h_bump_1_field__":0.0,
                                "__v_bump_1_field__":0.0,
                                "__h_bump_2_field__":0.0,
                                "__v_bump_2_field__":0.0,
                                "__h_bump_3_field__":0.0,
                                "__v_bump_3_field__":0.1,
                                "__h_bump_4_field__":0.0,
                                "__v_bump_4_field__":0.0,
                                "__h_bump_5_field__":0.0,
                                "__v_bump_5_field__":0.0,
                    },
                    "target_orbit":dict([
                        (0, [3897.68994, -17.78, 89.38, 3.927]), # WARNING this number is dubious
                        ]+[(i, co) for i in (1, 2,)]),
                    "position_tolerance":0.01,
                    "momentum_tolerance":0.001,
                    "fix_bumps":["__v_bump_-1_field__", "__h_bump_-1_field__",
                                 "__v_bump_-2_field__", "__h_bump_-2_field__",
                                 "__v_bump_3_field__", "__h_bump_3_field__",
                                 "__v_bump_4_field__", "__h_bump_4_field__",
                                 "__v_bump_5_field__", "__h_bump_5_field__",],
                }, { # recover closed orbit
                    "seed_fields":{},
                    "target_orbit":dict([ #maps station:target orbit
                        ]+[(i, co) for i in (1, 2, 4)]),
                    "position_tolerance":0.01,
                    "momentum_tolerance":0.001,
                    "fix_bumps":["__v_bump_1_field__", "__h_bump_1_field__",
                                 "__v_bump_2_field__", "__h_bump_2_field__",
                                 "__v_bump_3_field__", "__h_bump_3_field__",
                                 "__v_bump_-4_field__", "__h_bump_-4_field__",
                                 "__v_bump_-5_field__", "__h_bump_-5_field__",],
                }#
            ],
            "target_fields":{},
            "seed_errors":[1e-3]*10,
            "ref_probe_files":["FOILPROBE.h5", "RINGPROBE*.h5"], # sorted alphanumerically
            "run_dir":"tmp/find_bump/",
            "energy":3.0,
            "min_time_delta":0., # minimum time between probes
            "target_n_hits":3,
            "penalty_factor":1e9, # penalty = p_f^(number of missed stations)
            "algorithm":"migrad",
        }
        probes = [("FOILPROBE.h5", 0)]
        probes += [("RINGPROBE"+str(i).rjust(2, "0")+".h5", i) for i in range(1, 11)]
        probes += [("RINGHALFPROBE"+str(i).rjust(2, "0")+".h5", i+10) for i in range(1, 11)]
        probes = dict(probes)
        print("PROBES", probes)
        self.build_bump_surrogate_model = {
            "track_map":True,
            "build_map":True,
            "n_h_bumps":5,
            "n_v_bumps":5,
            "output_file":"build_bump_surrogate_model",
            "subs_overrides":{
                "__n_turns__":0.7,
                "__do_magnet_field_maps__":False,
                "__max_x_power__":10,
                "__do_bump__":True,
            },
            "closed_orbit":co,
            "ref_probe_files":probes, # sorted alphanumerically
            "energy":3.0, # sorted alphanumerically
            "run_dir":"tmp/build_bump_model/",
            "bump_probe_station":0,
            "fit_order":1,
            "smoothing_order":1,
            "staged_optimisation":[{
                    "name":"Upstream of foil",
                    "n_per_dimension":5, # number of field steps per dimension
                    "inverse_n_per_dimension":1, # inversion number of position steps per dimension
                    "field_scale":0.25, # fields extend +/- x [T]
                    "extrapolation_factor":101,
                    "seed_fields":{
                        "__h_bump_1_field__":0.0,
                        "__v_bump_1_field__":0.0,
                        "__h_bump_2_field__":0.0,
                        "__v_bump_2_field__":0.0,
                        "__h_bump_3_field__":0.0,
                        "__v_bump_3_field__":0.10,
                        "__h_bump_4_field__":0.0,
                        "__v_bump_4_field__":0.0,
                        "__h_bump_5_field__":0.0,
                        "__v_bump_5_field__":0.0,
                    },
                    "fix_bumps":[
                        "__v_bump_-1_field__", "__h_bump_-1_field__",
                        "__v_bump_-2_field__", "__h_bump_-2_field__",
                        "__v_bump_3_field__", "__h_bump_3_field__",
                        "__v_bump_4_field__", "__h_bump_4_field__",
                        "__v_bump_5_field__", "__h_bump_5_field__",
                    ],
                    "fore_tracking":True,
                    "phi_start":0.,
                },{
                    "name":"Downstream of foil",
                    "n_per_dimension":5, # number of field steps per dimension
                    "inverse_n_per_dimension":1, # inversion number of position steps per dimension
                    "field_scale":0.2, # fields extend +/- x [T]
                    "extrapolation_factor":101, # factor by which we are allowed to extrapolate fields
                    "seed_fields":{
                                "__h_bump_1_field__":0.0,
                                "__v_bump_1_field__":0.0,
                                "__h_bump_2_field__":0.0,
                                "__v_bump_2_field__":0.0,
                                "__h_bump_3_field__":0.0,
                                "__v_bump_3_field__":0.10,
                                "__h_bump_4_field__":0.0,
                                "__v_bump_4_field__":0.0,
                                "__h_bump_5_field__":0.0,
                                "__v_bump_5_field__":0.0,
                    },
                    "fix_bumps":["__v_bump_1_field__", "__h_bump_1_field__",
                                 "__v_bump_2_field__", "__h_bump_2_field__",
                                 "__v_bump_3_field__", "__h_bump_3_field__",
                                 "__v_bump_-4_field__", "__h_bump_-4_field__",
                                 "__v_bump_-5_field__", "__h_bump_-5_field__",],
                     "fore_tracking":False,
                     "phi_start":-4.,
                }],
        }
        self.track_bump = {
            "input_file":"find_bump_parameters.tmp",
            "injection_orbit":0, # reference to item from bump_list
            "subs_overrides":{
                "__n_turns__":1.9,
                "__no_field_maps__":"",
            },
            "bump_list":None, # list of bumps which we will track (default to find_bump_parameters)
            "n_turns_list":None, # number of turns for forward tracking (default 1)
            "foil_probe_files":["FOILPROBE.h5"], # PROBE where the beam is injected
            "ramp_probe_files":["RINGPROBE06.h5"], # PROBE where the magnet is ramped
            "ramp_probe_phi":5, # cell number where we start the beam following an injection
            "run_dir":"tmp/track_bump/",
            "energy":3.0,
            "bump_probe_station":0,
        }

        self.substitution_list = [get_baseline_substitution(offset)]
        
        self.run_control = {
            "find_closed_orbits_4d":False,
            "find_tune":False,
            "find_da":True,
            "find_bump_parameters":False,
            "build_bump_surrogate_model":False,
            "track_bump":False,
            "clean_output_dir":False,
            "output_dir":os.path.join(os.getcwd(), "output/double_triplet_baseline/offset_scan/"+str(offset).ljust(4, "0")+"/"),
            "root_verbose":6000,
            "faint_text":'\033[38;5;243m',
            "default_text":'\033[0m'
        }

        self.tracking = {
            "mpi_exe":None, #os.path.expandvars("${OPAL_ROOT_DIR}/external/install/bin/mpirun"),
            "beam_file_out":"disttest.dat",
            "n_cores":4,
            "links_folder":["VerticalSectorFFA",], # link relative to lattice/VerticalFFA.in
            "lattice_file":os.path.join(os.getcwd(), "lattice/VerticalDoubleTripletFFA.in"),
            "lattice_file_out":"VerticalSectorFFA.tmp",
            "opal_path":os.path.expandvars("${OPAL_EXE_PATH}/opal"),
            "tracking_log":"log",
            "flags":[],
            "step_size":1.,
            "ignore_events":[],
            "pdg_pid":2212,
            "clear_files":"*.h5",
            "verbose":0,
            "file_format":"hdf5",
            "analysis_coordinate_system":"azimuthal",
            "dt_tolerance":1., # ns
            "py_tracking":{
                "derivative_function":"u_derivative",
                "atol":tol,
                "rtol":tol,
                "verbose":False,
            }
        }
