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
        "__d_offset__":-0.154, # m
        "__f_offset__":0.154, # m
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
        "__magnet_width__":2.0, # m
        "__bb_length__":8.0, # m
        # field maps
        "__do_magnet_field_maps__":True,
        "__cartesian_x_min__":-5., # m
        "__cartesian_dx__":0.01, # m (1001 steps)
        "__cartesian_y_min__":-5., # m
        "__cartesian_dy__":0.01, # m (1001 steps)
    }
    return baseline

# Status 19/07/2019:
# * Tracking finds closed orbit using Shinji lattice
# * DA is tiny (0) and TM declares phase advance 0.5, 0.0
# Ideas
# * Try tracking backwards (with reversed polarities)
# * Try searching for good tune; first by tweaking end length, then by tweaking DF ratio

# Finder is failing; does not seem to be improved by adjusting step size. Check 
# the lattice? Try constraining fringe field and check probe field is 0 (in which
# case lattice should be perfect)?

class Config(object):
    def __init__(self):
        self.find_closed_orbits = { 
            "seed":[[3957.3378704328647, 22.748878261530315, -3.106790633028175, 1.0895278377208744]],
            # I tried with deltas [5., 0.05, 5., 0.05] and the TM determinant
            # converged on ~ 0.85 as a function of step size
            "deltas":[1., 0.01, 1., 0.01],
            "adapt_deltas":False,
            "output_file":"closed_orbits_cache",
            "subs_overrides":{"__n_turns__":0.21, "__do_magnet_field_maps__":"False"},
            "final_subs_overrides":{"__n_turns__":10.1, "__do_magnet_field_maps__":"True"},
            "root_batch":0,
            "max_iterations":5, 
            "tolerance":0.1,
            "ds_cell":2,
            "do_plot_orbit":False,
            "run_dir":"tmp2/",
            "probe_files":"RINGPROBE*.loss",
            "overwrite":True,
            "orbit_file":"VerticalFFAGMagnet-trackOrbit.dat",
        }
        self.find_tune = {
            "run_dir":"tmp/find_tune/",
            "probe_files":"RINGPROBE*.loss",
            "subs_overrides":{"__n_turns__":5.1, "__do_magnet_field_maps__":"True"},
            "root_batch":0,
            "delta_1":1.,
            "delta_2":1.,
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
            "opal_path":os.path.expandvars("${OPAL_EXE_PATH}/opal"),
            "tracking_log":"log",
            "step_size":1.,
            "pdg_pid":2212,
            "clear_files":"*.loss",
            "verbose":False,
        }


def get_baseline_substitution_rogers():
    baseline = {
        # ring
        "__n_cells__":10,
        "__radius__":3.995285, #m 
        "__drift_1__":0.55,
        "__drift_2__":0.55,
        # beam
        "__energy__":3., # MeV
        # tracking
        "__step_size__":0.001, # m
        "__n_turns__":0.2,
        # main magnets
        "__bf__":-0.43, # T
        "__bd__":0.2, #+0.20, # T
        "__d_length__":0.4, # m
        "__f_length__":0.4, # m
        "__d_offset__":-0.15, # m
        "__f_offset__":0.15, # m
        "__d_end_length__":0.125,  # m
        "__f_end_length__":0.125, # m
        "__m_index__":1.7, # m^-1
        "__max_x_power__":6,
        "__neg_extent__":0.2, # m
        "__pos_extent__":0.8, # m
        "__magnet_separation__":-0.1, # m
        "__f_tilt_angle__":-8.0,
        "__d_tilt_angle__":+4.0,
        # NOTE BB length 0.125 m -> width 1.0 m; length 2.0 m; double for safety
        "__magnet_width__":0.5, # m
        "__bb_length__":4.0, # m
        # field maps
        "__do_magnet_field_maps__":True,
        "__cartesian_x_min__":-5., # m
        "__cartesian_dx__":0.01, # m (1001 steps)
        "__cartesian_y_min__":-5., # m
        "__cartesian_dy__":0.01, # m (1001 steps)
    }
    return baseline

