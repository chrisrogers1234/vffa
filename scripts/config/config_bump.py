import math
import os

def get_baseline_substitution():
    baseline = {
        # ring
        "__n_cells__":10,
        "__radius__":3.995285, #m
        "__drift_1__":0.75,
        "__drift_2__":0.75,
        # beam
        "__energy__":3., # MeV
        # tracking
        "__step_size__":0.001, # m
        "__n_turns__":0.2,
        # main magnets
        "__bf__":-0.25/0.44, # T
        "__bd__":0.25, #+0.20, # T
        "__d_length__":0.2, # m
        "__f_length__":0.2, # m
        "__d_offset__":-0.154, # m # 0.154
        "__f_offset__":0.154, # m
        "__d_end_length__":0.125,  # m
        "__f_end_length__":0.125, # m
        "__m_index__":1.4, # m^-1 # 1.58
        "__max_x_power__":6,
        "__neg_extent__":2.0, # m
        "__pos_extent__":4.0, # m
        "__magnet_separation__":0.1, # m
        "__f_tilt_angle__":-16.0,
        "__d_tilt_angle__":+8.0,
        # NOTE BB length 0.125 m -> width 1.0 m; length 2.0 m; double for safety
        "__magnet_width__":0.5, # m
        "__bb_length__":8.0, # m
        # field maps
        "__do_magnet_field_maps__":True,
        "__cartesian_x_min__":-5., # m
        "__cartesian_dx__":0.01, # m (1001 steps)
        "__cartesian_y_min__":-5., # m
        "__cartesian_dy__":0.01, # m (1001 steps)
        # bump magnets
        "__do_bump__":True,
        "__bump_length__":0.1,
        "__bump_fringe__":0.1,
        "__bump_width__":2.0,
        "__bump_offset__":0.0, # radial offset, m
        "__bump_angle__":41.0, # angular tilt of bump cavity
        "__phi_foil_probe__":3.2,
        "__phi_bump_1__":+0.35,
        "__phi_bump_2__":-0.65,
        "__phi_bump_3__":-1.3,
        "__phi_bump_4__":-1.65,
        "__h_bump_1_field__":10.0,
        "__h_bump_2_field__":10.0,
        "__h_bump_3_field__":10.0,
        "__h_bump_4_field__":10.0,
        "__v_bump_1_field__":10.0,
        "__v_bump_2_field__":10.0,
        "__v_bump_3_field__":10.0,
        "__v_bump_4_field__":10.0,
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
            "seed":[[3932.3656651884976, -26.152597049937057, 28.78211467244023, -0.3868433565977967]],
            "deltas":[0.5, 0.5, 0.5, 0.5],
            "adapt_deltas":False,
            "output_file":"closed_orbits_cache",
            "subs_overrides":{"__n_turns__":0.21, "__do_magnet_field_maps__":"False"},
            "final_subs_overrides":{"__n_turns__":0.21, "__do_magnet_field_maps__":"True"},
            "root_batch":0,
            "max_iterations":0, 
            "tolerance":0.1,
            "ds_cell":2,
            "do_plot_orbit":False,
            "run_dir":"tmp/find_closed_orbits",
            "probe_files":"RINGPROBE*.loss",
            "overwrite":True,
            "orbit_file":"VerticalFFAGMagnet-trackOrbit.dat",
        }
        self.find_tune = {
            "run_dir":"tmp/find_tune/",
            "probe_files":"RINGPROBE*.loss",
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
            "probe_files":"RINGPROBE01.loss",
            "subs_overrides":{"__n_turns__":50.1, "__no_field_maps__":"// "},
            "get_output_file":"get_da",
            "scan_output_file":"scan_da",
            "row_list":None,
            "scan_x_list":[],
            "scan_y_list":[],
            "x_seed":1.,
            "y_seed":1.,
            "min_delta":1.0,
            "max_delta":500.,
            "required_n_hits":50,
            "max_iterations":10,
        }
        self.find_bump_parameters = {
            "n_h_bumps":4,
            "bump":[[0.0, 0.0, 0.0, 0.0 ] for i in range(1)],
            "output_file":"find_bump_parameters",
            "closed_orbit":[3932.3656651884976, -26.152597049937057, 28.78211467244023, -0.3868433565977967],
            "magnet_min_field":-1.0,
            "magnet_max_field":+1.0,
            "max_iterations":500,
            "position_tolerance":0.1,
            "momentum_tolerance":100.,
            "subs_overrides":{
                "__n_turns__":1.1,
                "__do_magnet_field_maps__":False,
            },
            "final_subs_overrides":{
                "__n_turns__":1.1,
                "__do_magnet_field_maps__":True,
            },
            "fix_bumps":["__v_bump_1__"],
            "seed_field":[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            "seed_errors":[0.1, 0.1, 0.1, 0.1],
            "foil_probe_phi":3.2,
            "bump_probe_station":0,
            "ignore_stations":[3, 4, 5, 6],
            "ref_probe_files":["RINGPROBE*.loss"],
            "run_dir":"tmp/find_bump/",
            "energy":3.0,
        }
        
        self.substitution_list = [get_baseline_substitution()]
        
        self.run_control = {
            "find_closed_orbits_4d":True,
            "find_tune":False,
            "find_da":False,
            "clean_output_dir":False,
            "output_dir":os.path.join(os.getcwd(), "output/baseline"),
            "root_verbose":6000,
        }

        self.tracking = {
            "mpi_exe":None, #os.path.expandvars("${OPAL_ROOT_DIR}/external/install/bin/mpirun"),
            "beam_file_out":"disttest.dat",
            "n_cores":4,
            "links_folder":["VerticalFFA",], # link relative to lattice/VerticalFFA.in
            "lattice_file":os.path.join(os.getcwd(), "lattice/VerticalFFA.in"),
            "lattice_file_out":"VerticalFFAGMagnet.tmp",
            "opal_path":os.path.expandvars("/home/vol08/scarf148/data/OPAL/opal_src_test/src/opal"), #"${OPAL_EXE_PATH}/opal"),
            "tracking_log":"log",
            "step_size":1.,
            "pdg_pid":2212,
            "clear_files":"*.loss",
            "verbose":False,
        }

