import math

import matplotlib
import numpy
import platypus
import json
from platypus import Problem
from toy_model import ToyModel

class ToyModelOptimisation(object):
    def __init__(self, config):
        self.model = ToyModel()
        self.model.do_plots = False
        #self.model.number_pulses = 50
        self.model.verbose = -1
        self.model.parse_args()
        if config != None:
            self.model.do_config(config)
        self.seed = [0, 0, 0]
        self.sigma = 1.
        self.n_trials = 1000
        self.final_event_multiplier = 10 # final run gets n times more per pulse
        self.iteration = 0
        self.algorithm = self.setup_optimisation_cmaes()
        self.iterations = []
        self.config = config

    def optimise(self):
        self.iteration = 0
        self.algorithm.run(self.n_trials)
        scores = sorted(self.iterations, key = lambda x: x[0][0])
        best = scores[0]
        print("Running best inputs", best[1], "with score", best[0])
        self.model.do_plots = True
        self.model.number_per_pulse *= self.final_event_multiplier
        self.run_one(best[1])
        self.finalise()
        #matplotlib.pyplot.show(block=False)

    def finalise(self):
        fout = open(self.model.output_dir+"/config.json", "w")
        print(json.dumps(self.config), file=fout)

    def correlated_default_bumps(self, angle_1, angle_2, target_emittance, n_pulses):
        self.model.default_bumps = [[angle_1, i*target_emittance/(n_pulses-1), # (n_pulses-1-i)*1.3/(n_pulses-1),
                                     angle_2, i*target_emittance/(n_pulses-1)]
                                     for i in range(n_pulses)]

    def anticorrelated_default_bumps(self, angle_1, angle_2, target_emittance, n_pulses):
        self.model.default_bumps = [[angle_1, (n_pulses-1-i)*target_emittance/(n_pulses-1),
                                     angle_2, i*target_emittance/(n_pulses-1)]
                                     for i in range(n_pulses)]

    def run_one(self, arguments):
        self.iteration += 1
        print("Run", self.iteration, "with args", arguments, end='')
        n_pulses = self.model.number_pulses
        angle_1 = arguments[0]
        angle_2 = arguments[1]
        foil_angle = math.degrees(arguments[2])
        self.anticorrelated_default_bumps(angle_1, angle_2, 1.3, n_pulses)
        self.model.foil_angle = foil_angle   
        self.model.initialise()
        self.model.run_toy_model()
        self.model.finalise()
        mean_foil_hits = numpy.mean(self.model.foil_hits)
        print("    ... Mean foil hits", mean_foil_hits)
        score = [mean_foil_hits] 
        self.iterations.append([score, arguments])
        return score

    def setup_optimisation_cmaes(self):
        n_parameters = 3
        n_variables = 1
        problem = platypus.Problem(n_parameters, n_variables)
        problem.types[:] = platypus.Real(-math.pi/2., math.pi/2.)
        problem.directions[:] = Problem.MINIMIZE
        problem.function = self.run_one
        algorithm = platypus.CMAES(problem)
        algorithm.initial_search_point = self.seed
        algorithm.sigma = self.sigma # note this is a single value
        return algorithm

def do_config(name, dp_over_p, col_density, emit):
    directory = "output/triplet_baseline/optimiser/"+name+"__"+\
                "dp_over_p_"+str(format(dp_over_p, "4.4g"))+"__"+\
                "col_dens_"+str(format(col_density, "4.4g"))+"__"+\
                "emit_"+str(format(emit, "4.4g")+"/")
    return {"dp_over_p":dp_over_p,
            "foil_column_density":col_density,
            "pulse_emittance":emit,
            "output_dir":directory,
            "number_per_pulse":10}

def main():
    name = "anticorrelated"
    config_list = [do_config(name, 1e-9, 1e-9, 1e-9),
                   do_config(name, 1e-9, 1e-9, 0.026),
                   do_config(name, 0.0013, 1e-9, 0.026),
                   do_config(name, 0.0013, 5e-6, 0.026),
                   do_config(name, 0.0013, 2e-5, 0.026)]
    for config in config_list:
        optimisation = ToyModelOptimisation(config)
        optimisation.optimise()

    matplotlib.pyplot.show(block=False)
    input("Press <CR> to finish")

if __name__ == "__main__":
    main()
