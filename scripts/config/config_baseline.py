import math
import os

def get_baseline_substitution():
    baseline = {
        # beam
        "__energy__":3., # MeV
        "__beam_theta__":0.,
        "__beam_charge__":1.,
        "__beam_phi_init__":0.0,
        "__beamfile__":'disttest.dat',
        "__current__":1.6e-19,
        "__n_events__":1,
        "__lattice_phi_offset__":0.0,
        # tracking
        "__step_size__":1., # mm
        "__n_turns__":2.1,
        "__solver__":"None",
        "__mx__":5,
        "__my__":5,
        "__mt__":5,
        "__test_probe__":0.0,
        # main magnets
        "__b_mean__":-3.*(1.*4.8-0.36*2.4)/24.,
        "__df_ratio__":-0.36,
        "__d_length__":0.1,
        "__f_length__":0.2,
        "__d_end_length__":(1.0/math.cos(math.radians(41)))/1.0/2.4,  #*0.9/2.,
        "__f_end_length__":(1.0/math.cos(math.radians(41)))/3.65/4.8,  #*0.5/2.,
        "__tan_delta__":-math.tan(math.radians(41)),
        "__field_index__":7.1,
        "__max_y_power__":2,
        "__neg_extent__":0.2,
        "__pos_extent__":0.8,
        # rf
        "__cavity__":"no_cavity",
        "__rf_voltage__":0.,
        "__rf_phase__":1.90862003999-math.radians(20), # [rad]
        "__rf_freq_0__":0.93656519779, # [MHz]
        "__rf_freq_1__":0.,
        "__rf_freq_2__":0.,
        "__rf_freq_3__":0.,
        "__no_field_maps__":"// ", # set to "" to enable field maps; set to "// " to comment them
        "__cavity_angle__":41.,
        # bump magnets
        "__bump_length__":0.1,
        "__bump_fringe__":0.1,
        "__bump_width__":1.0,
        "__bump_offset__":0.0, # radial offset, m
        "__bump_angle__":41.0, # angular tilt of bump cavity
        "__phi_foil_probe__":3.2,
        "__phi_bump_1__":+0.35,
        "__phi_bump_2__":-0.65,
        "__phi_bump_3__":-1.3,
        "__phi_bump_4__":-1.65,
        "__bump_field_1__":0.0, #-0.153622540713,
        "__bump_field_2__":0.0, #+0.0,
        "__bump_field_3__":0.0, #+0.157288331088,
        "__bump_field_4__":0.0, #-0.20832864219,
        "__septum_field__":0.0, #-0.20832864219,
    }
    return baseline

class Config(object):
    def __init__(self):
        self.find_closed_orbits = {
            "seed":[[4000.0, -8.47136893],],
            "output_file":"find_closed_orbit",
            "subs_overrides":{"__n_turns__":1.1, "__no_field_maps__":""},
            "root_batch":0,
            "max_iterations":10,
            "do_plot_orbit":False,
            "run_dir":"tmp/find_closed_orbits/",
            "probe_files":"RINGPROBE*.loss",
        }
        self.find_tune = {
            "run_dir":"tmp/find_tune/",
            "probe_files":"RINGPROBE01.loss",
            "subs_overrides":{"__n_turns__":50.1, "__no_field_maps__":"// ", "__step_size__":0.1},
            "root_batch":0,
            "delta_x":2.,
            "delta_y":1.,
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
            "find_closed_orbits":True,
            "find_tune":False,
            "find_da":False,
            "clean_output_dir":False,
            "output_dir":os.path.join(os.getcwd(), "output/baseline"),
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
        }

