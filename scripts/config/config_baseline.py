import math
import os

def get_baseline_substitution():
    baseline = {
        # ring
        "__n_cells__":10,
        "__radius__":3.995285, #m
        "__drift_1__":0.75,
        "__drift_2__":0.75,
        "__n_particles__":3,
        # beam
        "__energy__":3., # MeV
        "__beam_phi_init__":0., # fraction of cell length
        # tracking
        "__step_size__":0.001, # m
        "__n_turns__":0.2,
        # main magnets
        "__bf__":-0.5, #-0.25/0.44, # T
        "__bd__":0.22, #0.25, #+0.20, # T
        "__d_length__":0.2, # m
        "__f_length__":0.2, # m
        "__d_offset__":-0.1625, # m # 0.154
        "__f_offset__":0.1625, # m
        "__d_end_length__":0.125,  # m
        "__f_end_length__":0.125, # m
        "__m_index__":1.58, # m^-1 # 1.58
        "__max_x_power__":6,
        "__neg_extent__":2.0, # m
        "__pos_extent__":4.0, # m
        "__magnet_separation__":0.1, # m
        "__f_tilt_angle__":-8.0,
        "__d_tilt_angle__":+4.0,
        # NOTE BB length 0.125 m -> width 1.0 m; length 2.0 m; double for safety
        "__magnet_width__":6.0, # m
        "__bb_length__":20.0, # m
        # field maps
        "__do_magnet_field_maps__":True,
        "__cartesian_x_min__":-5., # m
        "__cartesian_dx__":0.01, # m (1001 steps)
        "__cartesian_y_min__":-5., # m
        "__cartesian_dy__":0.01, # m (1001 steps)
        # bump magnets
        "__do_bump__":False,
        "__foil_probe_phi__":2.07,
        "__h_bump_1_field__":0.0,
        "__h_bump_2_field__":0.0,
        "__h_bump_3_field__":0.0,
        "__h_bump_4_field__":0.0,
        "__v_bump_1_field__":0.0,
        "__v_bump_2_field__":0.0,
        "__v_bump_3_field__":0.0,
        "__v_bump_4_field__":0.0,
        "__h_bump_1_phi__":1.03, #1.43,
        "__h_bump_2_phi__":2.03, #2.43,
        "__h_bump_3_phi__":3.03, #3.43,
        "__h_bump_4_phi__":4.03, #4.43,
        "__v_bump_1_phi__":1.03,
        "__v_bump_2_phi__":2.03,
        "__v_bump_3_phi__":3.03,
        "__v_bump_4_phi__":4.03,
        "__septum_field__":0.0, #-0.20832864219,
        "__septum_length__":0.15, # m
        "__septum_fringe__":0.1, # fraction of septum_length
        "__septum_width__":0.5, # m
        "__septum_phi__":0.5, # fraction of cell length
        "__septum_dr__":0.5, # m
    }
    return baseline

class Config(object):
    def __init__(self):
        self.find_closed_orbits = { 
            "seed":[[3946.1223258417176, -24.618972507786452, 136.66277886171747, -1.235505893912702]],
            "deltas":[0.1, 0.1, 0.1, 0.1],
            "adapt_deltas":False,
            "output_file":"closed_orbits_cache",
            "subs_overrides":{"__n_turns__":0.21, "__do_magnet_field_maps__":"False"},
            "final_subs_overrides":{"__n_turns__":0.51, "__do_magnet_field_maps__":"True"},
            "us_cell":1,
            "ds_cell":2,
            "root_batch":0,
            "max_iterations":5, 
            "tolerance":0.1,
            "ds_cell":2,
            "do_plot_orbit":False,
            "run_dir":"tmp/find_closed_orbits",
            "probe_files":"RINGPROBE*.h5",
            "overwrite":True,
            "orbit_file":"VerticalFFA-trackOrbit.dat",
        }
        self.find_tune = {
            "run_dir":"tmp/find_tune/",
            "probe_files":"RINGPROBE*.h5",
            "subs_overrides":{"__n_turns__":10.1, "__do_magnet_field_maps__":"True"},
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
            "subs_overrides":{"__n_turns__":120.1, "__do_magnet_field_maps__":"False", "__step_size__":0.001},
            "get_output_file":"get_da",
            "scan_output_file":"scan_da",
            "row_list":None,
            "scan_x_list":[],
            "scan_y_list":[],
            "x_seed":32.,
            "y_seed":32.,
            "min_delta":0.9,
            "max_delta":500.,
            "required_n_hits":100,
            "dt_tolerance":0.5, # fraction of closed orbit dt
            "max_iterations":15,
            "decoupled":True,
        }
        self.find_bump_parameters = {
            "n_h_bumps":4,
            "n_v_bumps":4,
            "bump":[[0.0, 0.0, -10.0*i, 0.0] for i in range(5, 6)],
            "amp":[],
            #"bump_reference":[0.0, 0.0, 0.0, 0.0], # relative to closed orbit
            "output_file":"find_bump_parameters",
            "closed_orbit":[3946.1223258417176, -24.618972507786452, 136.66277886171747, -1.235505893912702],
            "foil_closed_orbit":[3890.56766, -21.37, 135.68, 0.78877],
            "magnet_min_field":-0.2,
            "magnet_max_field":+0.2,
            "max_iterations":1000,
            "position_tolerance":0.01,
            "momentum_tolerance":1000.,
            "field_tolerance":1e-4,
            "amplitude_tolerance":1.,
            "tm_source":"closed_orbits_cache",
            "subs_overrides":{
                "__n_turns__":1.1,
                "__do_magnet_field_maps__":False,
                "__do_bump__":True,
            },
            "final_subs_overrides":{
                "__n_turns__":1.1,
                "__do_magnet_field_maps__":True,
                "__do_bump__":True,
            },
            "fix_bumps":["__h_bump_2_field__"],
            "seed_field":[b*1. for b in
                          [-0.054195158303617735, 
                           0.0, 
                           -0.06153626386584041, 
                           0.006769920874463653, 
                           -0.012714180112733564, 
                           0.09356229195237886, 
                           -0.04667873523185467, 
                           0.09447524861712192]],
            "target_fields":{"__v_bump_2_field__":0.1},
            "seed_errors":[1e-3]*8,
            "bump_probe_station":0,
            "amp_probe_station":0,
            "ignore_stations":[3, 4, 5, 6],
            "ref_probe_files":["FOILPROBE.h5", "RINGPROBE*.h5"],
            "run_dir":"tmp/find_bump/",
            "energy":3.0,
            "min_time_delta":1., # minimum time between probes
            "target_n_hits":8,
            "penalty_factor":100, # penalty = p_f^(number of missed stations)
            "algorithm":"simplex",
        }
        self.track_bump = {
            "input_file":"find_bump_parameters.tmp",
            "injection_orbit":0, # reference to item from bump_list
            "subs_overrides":{
                "__n_turns__":1.2,
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
        
        self.run_control = {
            "find_closed_orbits_4d":False,
            "find_tune":False,
            "find_da":False,
            "find_bump_parameters":True,
            "track_bump":False,
            "clean_output_dir":False,
            "output_dir":os.path.join(os.getcwd(), "output/baseline_0.1625-3"),
            "root_verbose":6000,
        }

        self.tracking = {
            "mpi_exe":None, #os.path.expandvars("${OPAL_ROOT_DIR}/external/install/bin/mpirun"),
            "beam_file_out":"disttest.dat",
            "n_cores":4,
            "links_folder":["VerticalFFA",], # link relative to lattice/VerticalFFA.in
            "lattice_file":os.path.join(os.getcwd(), "lattice/VerticalFFA.in"),
            "lattice_file_out":"VerticalFFA.tmp",
            "opal_path":os.path.expandvars("${OPAL_EXE_PATH}/opal"),
            "tracking_log":"log",
            "flags":[],
            "step_size":1.,
            "ignore_events":[0],
            "pdg_pid":2212,
            "clear_files":"*.h5",
            "verbose":0,
            "file_format":"hdf5",
            "analysis_coordinate_system":"azimuthal",
        }

