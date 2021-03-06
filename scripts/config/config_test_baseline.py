import math
import os

def r0(n_cells, half_cell_length):
    phi = 2.*math.pi/n_cells/4.
    radius = half_cell_length/2./math.sin(phi)
    print("R0", radius)
    return radius

# Machida TUNES 0.232180291582827 -0.183879862734979
# Rogers TUNES 0.2389296203830 0.17677569135392243

def get_baseline_substitution():
    R0 = r0(15, 3)
    baseline = {
        # ring
        "__n_cells__":15,
        "__radius__":R0, #m
        "__drift_1__":1.8,
        "__drift_2__":1.8,
        # beam
        "__energy__":400., # MeV
        # tracking
        "__step_size__":0.001, # m
        "__n_turns__":0.2,
        # main magnets
        "__bf__":-3.0*1.64, # T
        "__bd__":+1.0*1.64, # T
        "__d_length__":1.2, # m
        "__f_length__":1.2, # m
        "__d_offset__":-0.347745, # m
        "__f_offset__":0.347745, # m
        "__d_end_length__":0.3,  # m
        "__f_end_length__":0.3, # m
        "__m_index__":0.87721, # m^-1
        "__max_x_power__":10,
        "__neg_extent__":1.0, # m
        "__pos_extent__":3.0, # m
        # NOTE BB length 0.125 m -> width 1.0 m; length 2.0 m; double for safety
        "__magnet_width__":2.0, # m
        "__bb_length__":8.0, # m
        # field maps
        "__do_magnet_field_maps__":True,
        "__cartesian_x_min__":-20.0, # m
        "__cartesian_dx__":40/1000., # m (1001 steps)
        "__cartesian_y_min__":-20.0, # m
        "__cartesian_dy__":40/1000., # m (1001 steps)
        #"__cartesian_x_min__":-6.0, # m
        #"__cartesian_dx__":6/1000., # m (1001 steps)
        #"__cartesian_y_min__":R0-3.0, # m
        #"__cartesian_dy__":6/1000., # m (1001 steps)
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

class Config(object):
    def __init__(self):
        self.find_closed_orbits = {
            "seed":[[14222.486846033247, 235.56577518528238, 67.33679096440049, 6.319732909645213]],
                #[[14219.523614290156, -0.2573059678153996*1000, 72.17318674617552, -0.0074218697584373106*1000]], #  631.531
            "deltas":[0.01, 0.01, 0.01, 0.01],
            "adapt_deltas":False,
            "output_file":"closed_orbits_cache",
            "subs_overrides":{"__n_turns__":0.11, "__do_magnet_field_maps__":"False"},
            "final_subs_overrides":{"__n_turns__":0.21, "__do_magnet_field_maps__":"True"},
            "root_batch":0,
            "max_iterations":5,
            "tolerance":0.01,
            "ds_cell":1,
            "do_plot_orbit":False,
            "run_dir":"tmp/find_closed_orbits",
            "probe_files":"TESTPROBE*.loss",
            "overwrite":True,
            "orbit_file":"TestRing-trackOrbit.dat",
        }
        self.find_tune = {
            "run_dir":"tmp/find_tune/",
            "probe_files":"RINGPROBE*.loss",
            "subs_overrides":{"__n_turns__":10.1, "__do_magnet_field_maps__":"True"},
            "root_batch":0,
            "delta_1":5,
            "delta_2":5,
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
            "find_bump_parameters":False,
            "track_bump":False,
            "clean_output_dir":False,
            "output_dir":os.path.join(os.getcwd(), "output/test_baseline"),
            "root_verbose":6000,
        }

        self.tracking = {
            "mpi_exe":None, #os.path.expandvars("${OPAL_ROOT_DIR}/external/install/bin/mpirun"),
            "beam_file_out":"disttest.dat",
            "n_cores":4,
            "links_folder":["VerticalFFA",], # link relative to lattice/VerticalFFA.in
            "lattice_file":os.path.join(os.getcwd(), "lattice/JanFFA.in"),
            "lattice_file_out":"VerticalFFAMagnet.tmp",
            "opal_path":os.path.expandvars("${OPAL_BUILD_PATH}/tests/opal_unit_tests"),
            "flags":["--gtest_filter=VerticalFFAMagnetTrackingTest.*"],
            "tracking_log":"log",
            "step_size":1.,
            "pdg_pid":2212,
            "clear_files":None,
            "verbose":False,
            "file_format":"ascii",
        }


