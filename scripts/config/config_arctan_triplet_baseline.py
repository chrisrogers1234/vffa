import json
import math
import os

def get_baseline_substitution():
    MHz = 1e3
    injection_tof = 1151.534795
    field_factor = 0.22*math.exp(-1.31*0.116713)
    baseline = {
        # ring
        "__n_cells__":10,
        "__radius__":4.530495168499706, #m
        "__half_cell_length__":1.4,
        "__n_particles__":3,
        "__drift_radius__":4.357, # approx position where the beam passes through drift centre
        "__injection_tof__":injection_tof, # time of flight for a single turn
        # beam
        "__energy__":3., # MeV
        "__beam_phi_init__":0., # fraction of cell length
        "__beam_charge__":1., # multiples of positron charge
        # tracking
        "__step_size__":0.01, # m
        "__n_turns__":0.2,
        "__hdf5__":True,
        "__stepper__":"RK-4",
        "__spt_frequency__":1,
        "__output_algorithm__":"RK4",
        "__output_plane_tolerance__":1e-6,
        # main magnets
        "__end_field_model__":"arctan",
        "__bf__":-2.0*field_factor, #-0.25/0.44, # T
        "__bf_offset__":0.0,
        "__bd__":2.3*field_factor, #0.25, #+0.20, # T
        "__bd_offset__":0.0,
        "__d_length__":0.5, # m
        "__f_length__":0.5, # m
        "__fd_gap__":0.03, # m # 0.08
        "__drift_length__":1.24, # m
        "__f_tilt_angle__":0.0,  # m
        "__m_index__":-1.31, # m^-1 # 1.28
        "__d_end_length__":0.15,  # m
        "__f_end_length__":0.15, # m
        "__max_x_power__":10,
        "__neg_extent__":5.0, # m
        "__pos_extent__":5.0, # m
        "__bb_azimuthal__":54.0, # azimuthal extent in degrees relative to the centre of the current cell
        # NOTE BB length 0.125 m -> width 1.0 m; length 2.0 m; double for safety
        "__magnet_width__":10.0, #2*math.sin(math.radians(36))*2.8, #6.0, # m
        "__bb_length__":2*(1.4+2.8*math.cos(math.radians(36))), # m
        # field maps
        "__do_magnet_field_maps__":False,
        "__cartesian_x_min__":-5., # m
        "__cartesian_dx__":0.01, # m (1001 steps)
        "__cartesian_y_min__":-5., # m
        "__cartesian_dy__":0.01, # m (1001 steps)
        # rf
        "__do_rf__":False,
        "__rf_phi__":5.0,
        "__rf_length__":0.2,
        "__rf_end_factor__":0.2,
        "__rf_fringe_order__":0,
        "__do_rf_fringe_field__":True,
        "__rf_time__":[-1e9, 30*injection_tof, 1e6, 1e9],
        "__rf_frequency__":[2.0/injection_tof*MHz, 2.0/injection_tof*MHz, 2.0/injection_tof*MHz, 2.0/injection_tof*MHz],
        "__rf_efield__":[0.02]*4,
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

        "__quad_1_phi__":0.0,
        "__quad_1_dphi__":0.0,
        "__quad_1_grad__":0.0,
        "__quad_1_scales__":[1.0, 1.0],
        "__quad_1_times__":[-1e9, 1e9],

        "__foil_probe_1_phi__":3.095,
        "__foil_probe_1_dphi__":0.25/2.+0.025*math.pi/10.0,
        "__foil_probe_2_phi__":2.905,
        "__foil_probe_2_dphi__":0.25/2.+0.025*math.pi/10.0,
        "__foil_phi__":0.01,
        "__foil_dphi__":0.0,
        "__foil_min_z__":-0.10,
        "__foil_max_z__":-0.15,
        "__foil_inner_radius__":4.35,
        "__foil_outer_radius__":4.40,
        "__foil_thickness__":10e-6,
        "__do_foil__":False,

        "__septum_field__":-0.0,
        "__septum_length__":0.1, # m
        "__septum_fringe__":0.1, # fraction of septum_length
        "__septum_width__":0.5, # m
        "__septum_phi__":0.5, # fraction of cell length
        "__septum_dr__":0.5, # m
    }
    return baseline

class Config(object):
    """
    Adjust the configuration to use the correct bounding box; make step size
    small in the CO finder
    with large step size DA in u was 
    Axis x Seed 1284.375 Number of cells passing time cut/total hits 27 / 27 in 13.439568758010864 [s]
    Tried DA with step size 0.001 in the DA finder; no improvement (killed because I got bored)
    """

    def __init__(self):
        delta = 0.01
        xp_ratio = 0.001
        tol = 1e-12
        seed_list = []
        for xi in range(3):
            for yi in range(6):
                seed_list.append([4300.0+xi*100.0, 0.0, yi*100.0, 0.0])
        self.ring_tof = 1151.534795

        self.find_closed_orbits = {
            "seed":[[4357.790553052024, 0.0, 0.0, 0.0]], #116.71346624236696, 0.0]], # [4348.98, 0, -103.172, 0]],# end length 0.15, bb width exactly covers next cell 2*2.8*sin(36)
            "deltas":[delta, delta*xp_ratio, delta, delta*xp_ratio],
            "adapt_deltas":False,
            "output_file":"closed_orbits_cache",
            "subs_overrides":{"__n_turns__":0.11, "__do_magnet_field_maps__":"False", "__do_bump__":True, "__step_size__":0.01},
            "final_subs_overrides":{"__n_turns__":4.01, "__do_magnet_field_maps__":"False",
                    "__do_bump__":True, "__do_rf__":False, "__step_size__":0.01},
            "plotting_subs":{},
            "us_cell":0,
            "ds_cell":1,
            "root_batch":0,
            "max_iterations":0,
            "tolerance":0.00001,
            "do_plot_orbit":False,
            "do_minuit":True,
            "minuit_weighting":{"x":1, "x'":0.0001, "y":1, "y'":0.0001},
            "minuit_tolerance":1e-10,
            "minuit_iterations":100,
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
            "probe_files":"RINGPROBE01.h5",
            "subs_overrides":{"__n_turns__":26.1, "__do_magnet_field_maps__":"False", "__step_size__":0.01},
            "get_output_file":"get_da_2",
            "scan_output_file":"scan_da",
            "row_list":None,
            "scan_x_list":[],
            "scan_y_list":[],
            "x_seed":2400.0,
            "y_seed":None,
            "min_delta":0.9,
            "max_delta":1000.,
            "required_n_hits":25,
            "dt_tolerance":0.5, # fraction of closed orbit dt
            "max_iterations":10,
            "decoupled":True,
            "save_dir":"find_da",
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

        self.substitution_list = [get_baseline_substitution()]

        T0 = self.ring_tof
        self.track_beam = {
            "run_dir":"tmp/track_beam/",
            "save_dir":"track_beam_rf_on",
            "print_events":[i for i in range(1)],
            "variables":["x", "x'", "y", "y'", "t", "energy"],
            "settings":[{
                "name":"grid",
                "direction":"forwards",
                "probe_files":"RINGPROBE01.h5",          
                "beam":{
                    "type":"grid",
                    "energy":3.0,
                    "start":[4357.646683446333, 0.0, 0.0, 0.0, T0*0.5-0.125, 941.272],
                    "stop":[4357.646683446333, 0.0, 0.0, 0.0, T0*0.5+0.125, 941.272],
                    "nsteps":[1, 1, 1, 1, 40+1, 1],
                },
                "subs_overrides":{
                    "__n_turns__":20.1,
                    "__hdf5__":"True",
                    "__do_magnet_field_maps__":"False",
                    "__do_bump__":"False",
                    "__do_rf__":True,
                    "__step_size__":0.01
                },
            },],
        }
        max_amp = 0.02
        self.track_beam_null = {
            "run_dir":"tmp/track_beam/",
            "save_dir":"track_beam_amplitude_3",
            "print_events":[0, 1, -1],
            "settings":[{
                "name":"forwards",
                "direction":"forwards",
                "probe_files":"RINGPROBE01.h5",          
                "beam":{
                    "type":"beam_gen",
                    "closed_orbit_file":"closed_orbits_cache",
                    "eigen_emittances":[[0, 0]]*2+[[0.0020, 0], [0.02, 0]],
                    "n_per_dimension":4,
                    "variables":["x","x'","y","y'"],
                    "amplitude_dist":"grid", #"grid", # 
                    "phase_dist":"grid", #"grid", # 
                    "max_amplitude_4d":0.1, # amplitude_dist != grid
                    "energy":3.0,
                },
                "subs_overrides":{
                    "__n_turns__":1.1,
                    "__hdf5__":"True",
                    "__do_magnet_field_maps__":"False",
                    "__do_bump__":"False",
                    "__step_size__":0.01
                },
            },],
        }
        
        self.run_control = {
            "find_closed_orbits_4d":True,
            "find_tune":False,
            "find_da":False,
            "find_bump_parameters":False,
            "build_bump_surrogate_model":False,
            "track_bump":False,
            "track_beam":False,
            "clean_output_dir":False,
            "output_dir":os.path.join(os.getcwd(), "output/arctan_baseline/baseline"),
            "ffa_version":"2021-09-15 arctan vffa",
            "injection_concept":"horizontal_dispersive",
            "injection_version":"v1a",
            "root_verbose":6000,
            "faint_text":'\033[38;5;243m',
            "default_text":'\033[0m',
            "random_seed":0,
        }

        self.tracking = {
            "mpi_exe":None, #os.path.expandvars("${OPAL_ROOT_DIR}/external/install/bin/mpirun"),
            "beam_file_out":"disttest.dat",
            "n_cores":4,
            "links_folder":["VerticalSectorFFA",], # link relative to lattice/VerticalFFA.in
            "lattice_file":os.path.join(os.getcwd(), "lattice/VerticalTripletFFA.in"),
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
            "station_dt_tolerance":-1., # ns, if +ve and two hits are close together, reallocate station
            "py_tracking":{
                "derivative_function":"u_derivative",
                "atol":tol,
                "rtol":tol,
                "verbose":True,
            }
        }


    def softcoded_bump_fields(self, file_name):
        fin = open(file_name)
        last = {"score":None}
        for line in fin.readlines():
            json_data = json.loads(line)
            if last["score"] == None or last["score"] > json_data["score"]:
                last = json_data
        overrides = last["bump_fields"]
        overrides["__do_bump__"] = True
        self.print_bump_fields(overrides)
        return overrides

    def zero_bump_fields(self, subs):
        for key in subs:
            if "bump" in key and "field" in key:
                subs[key] = 0.0
        subs["__do_bump__"] = True # Set to true to get foil probe...
        return subs

    def sort_list(self, file_list):
        file_list_tmp = []
        for name in file_list:
            r = name.split("find_bump_r_")[1]
            r = r.split("_theta")[0]
            r = float(r)
            theta = name.split("_theta_")[1]
            theta = theta.split("/find_bump")[0]
            theta = float(theta)
            file_list_tmp.append([(r, theta), name])
        file_list_tmp = sorted(file_list_tmp)
        file_list_tmp = [f for f in reversed(file_list_tmp)]
        return file_list_tmp

    def ramp_fields(self, subs, file_list, turn_list, ring_tof, will_step):
        """
        subs: a dictionary; bump field settings will be appended or overwritten
        file_list: list of files from which to pull bump settings. file_list will
              be sorted using self.sort_list
        turn_list: list of turn numbers at which bump settings will be applied
        ring_tof: ring time-of-flight. Used to calculate time for each turn
        will_step: set to true to step the bump magnet settings rather than
              smoothly varying them.
        """
        json_data_list = []
        for nominal_coordinates, file_name in file_list:
            fin = open(file_name)
            best = None
            for line in fin.readlines():
                json_line = json.loads(line)
                if not best or json_line["score"] < best["score"]:
                    best = json_line
            print("Adding beam at", nominal_coordinates, "from file", file_name, "with score", best["score"], "Tracking", [t for t in best["tracking"] if t[0] == 4] )
            json_data_list.append(best)

        overrides = {}
        for file_index, json_line in enumerate(json_data_list):
            for key, value in json_line["subs"].items():
                if "_bump_" not in key or "_field_" not in key:
                    continue
                scale_key = key.replace("field", "scales")
                time_key = key.replace("field", "times")
                azimuth_key = key.replace("field", "phi")
                azimuth = subs[azimuth_key]/subs["__n_cells__"]
                if scale_key not in overrides:
                    overrides[key] = 1.0 #scale from 1.0
                    overrides[scale_key] = [value]
                    overrides[time_key] = [-1e9]
                if will_step:
                    overrides[scale_key] += [value, value]
                    overrides[time_key] += [(file_index-0.25)*ring_tof+1e-9, (file_index+0.75)*ring_tof]
                else:
                    time = (turn_list[file_index]+azimuth)*ring_tof #there is a phase shift owing to the azimuthal angle of the particle
                    overrides[scale_key] += [value]
                    overrides[time_key] += [time]
        for key in overrides.keys():
            if "bump" in key and "scales" in key:
                overrides[key] += [0.0, 0.0]
            elif "bump" in key and "times" in key:
                overrides[key] += [overrides[key][-1]+ring_tof, 1e9]
        subs.update(overrides)
        return subs

    def print_bump_fields(self, subs):
        print_str = ""
        print_dict = {}
        for key, value in subs.items():
            if "_bump_" in key and "_field__" in key:
                print_str += "    "+key.replace("__", "").replace("_", " ")+"   "+format(value, "8.5g")+"\n"
                print_dict[key] = subs[key]
        print(json.dumps(print_dict, indent=2))
        print(print_str)

