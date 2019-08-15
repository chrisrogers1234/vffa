import math
import os

def r0(n_cells, half_cell_length):
    phi = 2.*math.pi/n_cells/4.
    radius = half_cell_length/2./math.sin(phi)
    print "R0", radius
    return radius


def get_baseline_substitution():
    R0 = r0(15, 3)
    baseline = {
        # ring
        "__n_cells__":15,
        "__radius__":R0, #m
        "__drift_1__":2.4,
        "__drift_2__":2.4,
        # beam
        "__energy__":400., # MeV
        # tracking
        "__step_size__":0.001, # m
        "__n_turns__":0.2,
        # main magnets
        "__bf__":-3.0, # T
        "__bd__":+1.0, # T
        "__d_length__":0.6, # m
        "__f_length__":0.6, # m
        "__d_offset__":-0.347745, # m
        "__f_offset__":0.347745, # m
        "__d_end_length__":0.3,  # m
        "__f_end_length__":0.3, # m
        "__m_index__":0.877262, # m^-1 # 1.58
        "__max_x_power__":6,
        "__neg_extent__":1.0, # m
        "__pos_extent__":3.0, # m
        # NOTE BB length 0.125 m -> width 1.0 m; length 2.0 m; double for safety
        "__magnet_width__":1.0, # m
        "__bb_length__":4.0, # m
        # field maps
        "__do_magnet_field_maps__":True,
        #"__cartesian_x_min__":-20.0, # m
        #"__cartesian_dx__":40/1000., # m (1001 steps)
        #"__cartesian_y_min__":-20.0, # m
        #"__cartesian_dy__":40/1000., # m (1001 steps)
        "__cartesian_x_min__":-6.0, # m
        "__cartesian_dx__":6/1000., # m (1001 steps)
        "__cartesian_y_min__":R0-3.0, # m
        "__cartesian_dy__":6/1000., # m (1001 steps)
    }
    return baseline


seeds = []
for ri in range(7):
    for rpi in range(3):
        for zi in range(3):
            r = 14000.+50.*ri
            rp = -0.1-0.1*rpi
            z = 1150.+50.*zi
            seeds.append([r, rp, z, 0.0])
#14250         -0          0          0

class Config(object):
    def __init__(self):
        self.find_closed_orbits = {
            "seed":[[14224., -0.227, 1484.49, -0.00749917]],
            # I tried with deltas [5., 0.05, 5., 0.05] and the TM determinant
            # converged on ~ 0.85 as a function of step size
            "deltas":[1., 0.01, 1., 0.01],
            "adapt_deltas":False,
            "output_file":"closed_orbits_cache",
            "subs_overrides":{"__n_turns__":0.21, "__do_magnet_field_maps__":"False"},
            "final_subs_overrides":{"__n_turns__":10.1, "__do_magnet_field_maps__":"True"},
            "root_batch":0,
            "max_iterations":0,  
            "tolerance":0.1,
            "ds_cell":2,
            "do_plot_orbit":False,
            "run_dir":"tmp/",
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
            "output_dir":os.path.join(os.getcwd(), "output/jan_baseline"),
            "root_verbose":6000,
        }

        self.tracking = {
            "mpi_exe":None, #os.path.expandvars("${OPAL_ROOT_DIR}/external/install/bin/mpirun"),
            "beam_file_out":"disttest.dat",
            "n_cores":4,
            "links_folder":["VerticalFFA",], # link relative to lattice/VerticalFFA.in
            "lattice_file":os.path.join(os.getcwd(), "lattice/JanFFA.in"),
            "lattice_file_out":"VerticalFFAGMagnet.tmp",
            "opal_path":os.path.expandvars("${OPAL_EXE_PATH}/opal"),
            "tracking_log":"log",
            "step_size":1.,
            "pdg_pid":2212,
            "clear_files":"*.loss",
            "verbose":False,
        }


