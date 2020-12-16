import math
import json
import matplotlib
from _toy_model import ToyModel

def toy_models_single_turn_injection(foil_n_turns_list, a_dir):
    for foil, n_turns in foil_n_turns_list:
        t_list = [n for n in range(n_turns)]
        output_dir = a_dir+"/n_"+str(n_turns)+"__foil_"+str(foil)+"/"
        config = {
            "max_turn":n_turns,
            "turn_bumper_index":[t_list, t_list],
            "default_bumps":[["coupled", 0.0, 0.0, -40.0*i/(n_turns-1), 0.0] for i in t_list],
            "plot_frequency":n_turns+1,
            "foil_column_density":foil,
            "output_dir":output_dir,
            "number_per_pulse":1000,
            "closed_orbit":"output/triplet_baseline/baseline/closed_orbits_cache",
            "injection_ellipse_algorithm":"from_twiss",
            "beta_x":2.0,
            "alpha_x":1.0,
            "beta_y":0.5,
            "alpha_y":-1.0,
        }
        yield config

def toy_models_correlated_painting(a_dir, angle_list):
    for angle_u, angle_v in angle_list:
        n_injection_turns = 10
        n_trajectory_turns = 10
        n_turns = n_injection_turns+n_trajectory_turns
        t_list = [n for n in range(n_turns)]
        output_dir = a_dir+"/test_u_"+str(math.degrees(angle_u))+"_v_"+str(math.degrees(angle_v))+"/"
        foil = 5e-6
        final_amplitude = 2.0
        bumps = [["action_angle", angle_u, i*final_amplitude/(n_turns-1),
                  angle_v, i*final_amplitude/(n_turns-1)] for i in range(n_turns)]
        config = {
            "max_turn":n_turns,
            "turn_bumper_index":[t_list, t_list],
            "default_bumps":bumps,
            "plot_frequency":n_turns+1,
            "foil_column_density":foil,
            "output_dir":output_dir,
            "number_pulses":n_injection_turns,
            "number_per_pulse":1000,
            "closed_orbit":"output/triplet_baseline/baseline/closed_orbits_cache",
            "injection_ellipse_algorithm":"from_twiss",
            "beta_x":2.0,
            "alpha_x":1.0,
            "beta_y":0.5,
            "alpha_y":-1.0,
        }
        yield config


def run_one(config):
    model = ToyModel()
    #model.verbose = -1
    model.parse_args()
    if config != None:
        model.do_config(config)
    model.initialise()
    model.run_toy_model()
    params = model.get_output_parameters()
    model.finalise()
    return params

def save_output(output, a_dir):
    for out_dict in output:
        for key, value in out_dict.items():
            print(key, format(value, "6.4g"), end="  ") 
        print()
    fout = open(a_dir+"/run_summary_list.json", "w")
    fout.write(json.dumps(output))

def main():
    a_dir = "output/triplet_baseline/test"
    output = []
    for i, config in enumerate(toy_models_single_turn_injection(a_dir, foil_n_turns_list)):
        out_dict = run_one(config)
        output.append(out_dict)
    save_output(output, a_dir)
    matplotlib.pyplot.show(block=False)
    input("Press <CR> to finish")

if __name__ == "__main__":
    main()

