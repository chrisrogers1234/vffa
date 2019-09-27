import json
import sys
import math
import copy
import os
import shutil
import time

import ROOT
import numpy
numpy.seterr('raise')
numpy.set_printoptions(linewidth=200)

from xboa.hit import Hit
import utils.utilities
import xboa.common


from PyOpal.polynomial_coefficient import PolynomialCoefficient
from PyOpal.polynomial_map import PolynomialMap
from utils.decoupled_transfer_matrix import DecoupledTransferMatrix

class ClosedOrbitFinder4D(object):
    def __init__(self, config):
        self.config = config
        self.config_co = config.find_closed_orbits
        self.energy = self.config.substitution_list[0]["__energy__"]
        self.run_dir = os.path.join(self.config.run_control["output_dir"],
                                    self.config_co["run_dir"])
        self.var_list = ["x", "px", "y", "py"]
        self.subs_tracking = {}
        self.output_list = []
        self.print_list = []
        self.get_tracking(True)

    def get_tracking(self, clear):
        self.here = os.getcwd()
        probes = self.config_co["probe_files"]
        if clear:
            utils.utilities.clear_dir(self.run_dir)
        self.tracking = utils.utilities.setup_tracking(self.config, probes, self.energy)

    def seed_to_hit(self, seed, t):
        hit = utils.utilities.reference(self.config, self.energy)
        for i, var in enumerate(self.var_list):
            hit[var] = seed[i]
        hit["t"] = t
        hit.mass_shell_condition("pz") # adjust pz so E^2 = p^2 + m^2
        return hit

    def track_one(self, seed, t, is_final):
        overrides = self.config_co["subs_overrides"]
        if is_final:
            overrides = self.config_co["final_subs_overrides"]
        hit = self.seed_to_hit(seed, t)
        #print 'energy', hit['kinetic_energy'], 'p', hit['p'], 'pvec', hit['px'], hit['py'], hit['pz']
        os.chdir(self.run_dir)
        self.subs_tracking = utils.utilities.do_lattice(self.config, self.subs, overrides)
        track = self.tracking.track_one(hit)
        os.chdir(self.here)
        return track

    def get_decoupled(self, tm):
        m = tm.get_coefficients_as_matrix()
        dim = len(m)
        for i, row in enumerate(m):
            row = row[1:5]
            m[i] = row
        DecoupledTransferMatrix.det_tolerance = 1e9
        decoupled = DecoupledTransferMatrix(m)
        return decoupled

    def get_co(self, tm):
        # \vec{x} = \vec{m} + \matrix{M} \vec{x}
        # Solve (\matrix{1} - \matrix{M}) \vec{x} = \vec{m}
        m = tm.get_coefficients_as_matrix()
        print("TM:")
        for i, row in enumerate(m):
            m[i] = row[0:5]
        print(self.str_matrix(m))
        try:
            dm = self.get_decoupled(tm)
            print("Determinant:  ", numpy.linalg.det(dm.m))
            print("Phase advance:", [dm.get_phase_advance(i)/math.pi/2. for i in range(2)])
        except Exception:
            sys.excepthook(*sys.exc_info())
            dm = DecoupledTransferMatrix(numpy.identity(4))
            print("Failed to calculate phase advance")
        dim = len(m)
        m_vec = numpy.array([m[i][0] for i in range(dim)])
        matrix = numpy.array([[-m[i][j] for j in range(1, dim+1)] for i in range(dim)])
        for i in range(dim):
            matrix[i, i] += 1.
        inv_matrix = numpy.linalg.inv(matrix)
        x_vec = numpy.dot(inv_matrix, m_vec.transpose())
        m_orig = numpy.array([[m[i][j] for j in range(1, dim+1)] for i in range(dim)])
        x_vec = [x_vec[i] for i in range(dim)]
        return x_vec, dm

    def str_matrix(self, matrix, fmt="10.6g"):
        output_string = ""
        for item in matrix:
            try:
               output_string += self.str_matrix(item, fmt)+"\n"
            except TypeError:
               output_string += format(item, fmt)+" "
        return output_string

    def get_index_1(self, dim):
        yield [1]*dim
        index = [0]*dim
        index[0] = -1
        while sum(index) < len(index)*2:
            # update index
            for i, ind in enumerate(index):
                if ind == 2:
                    index[i] = 0
                else:
                    index[i] += 1
                    break
            yield index

    def get_index_2(self, dim):
        index = [1]*dim
        yield index
        for i in range(4):
            index = [1]*dim
            index[i] = 0
            yield index
        for i in range(4):
            index = [1]*dim
            index[i] = 2
            yield index

    def get_values(self, dim): # makes a hypercube extending from 0 through 1
        for i in [0, 1]:
            if dim == 1:
                yield[i]
            else:
                for value in self.get_values(dim-1):
                    yield [i]+value

    def get_values_2(self, dim): # makes a hyperdiamond at +- 1
        value_list = [[0. for i in range(dim)]]
        for d in range(dim):
            for delta in [-1, 1]:
                value = [0. for i in range(dim)]
                value[d] = delta
                value_list.append(value)
        for value in value_list:
            yield value

    def fit_matrix(self, tm_list_in, tm_list_out):
        coefficients = []
        print("fit_matrix TM LIST")
        print("    ", tm_list_in[0])
        print("    ", tm_list_out[0])
        for i in range(50):
            index = PolynomialMap.index_by_power(i, 4)
            do_continue = False
            do_break = False
            sum_index = []
            for i, n in enumerate(index):
                if n > 1:
                    do_break = True
                sum_index += [i]*n 
            if len(sum_index) < 2:
                do_continue = True
            #print i, "**", do_continue, do_break, "**", index, "**", sum_index
            if do_break:
                break
            if do_continue:
                continue
            coefficients += [PolynomialCoefficient(sum_index, j, 0.0) for j in range(4)]
        tm = PolynomialMap.least_squares(tm_list_in, tm_list_out, 1, coefficients)
        return tm

    def get_tm(self, seeds, deltas):
        ds_cell = self.config_co["ds_cell"]
        dim = len(seeds)
        ref_track = None
        tm_list_in, tm_list_out = [], []
        for j, values in enumerate(self.get_values_2(dim)):
            values = [seeds[i]+values[i]*deltas[i] for i in range(dim)]
            a_track = self.track_one(values, 0., False)
            print("\r", "-/|\\"[j%4], j, end=' ')
            sys.stdout.flush()
            if len(a_track) < ds_cell+1:
                raise RuntimeError("Output had "+str(len(a_track))+\
                            " points which is less than ds_cell+1 "+str(ds_cell+1))
            vhit_1 = [a_track[0][var]-seeds[i] for i, var in enumerate(self.var_list)]
            vhit_2 = [a_track[ds_cell][var]-seeds[i] for i, var in enumerate(self.var_list)]
            tm_list_in.append(vhit_1)
            tm_list_out.append(vhit_2)
            if j == 0:
                ref_track = a_track
        ref_track_in = [tm_list_in[0][i]+x for i, x in enumerate(seeds)]
        ref_track_out = [tm_list_out[0][i]+x for i, x in enumerate(seeds)]
        return tm_list_in, tm_list_out, ref_track

    def check_stuck(self, errors, new_co, tm_list_in, tm_list_out, seeds):
        tolerance = self.config_co["tolerance"]
        indices = list(range(len(new_co)))
        ref_in = tm_list_in[0]
        ref_out = tm_list_in[1]
        if len(errors) > 1:
            delta_error = errors[-2] - errors[-1]
            print("Checking stuck; error", errors[-1], "delta", delta_error, "tolerance", tolerance)
            if delta_error < tolerance and errors[-1] > tolerance:
                printout = [seeds[i]+new_co[i] for i in indices]
                print("Detected stuck finder a; giving it a kick from\n   ", printout, end=' ')
                new_co = [(ref_in[i]+ref_out[i])/2. for i in indices]
                printout = [seeds[i]+new_co[i] for i in indices]
                print("to\n   ", printout)
                errors[-1] = errors[0]
                print("and resetting last error to", errors[-1])
        return new_co

    def get_error(self, delta):
        scale = self.config_co["deltas"]
        error = sum([(delta[i]/scale[i])**2 for i in range(len(delta))])
        error = error**0.5
        return error

    def print_ref_track(self, ref_track, seeds, tm):
        print("Reference track for seed", seeds)
        for i, hit in enumerate(ref_track):
            ref_list = [hit[var] for var in self.var_list]
            fmt = "10.6g"
            if not tm:
                fmt = "14.10g"
            print("Ref track:", self.str_matrix(ref_list, fmt), end=' ')
            if not tm:
                print()
                continue
            print("decoupled", end=' ')
            coupled = [hit[var]-seeds[j] for j, var in enumerate(self.var_list)]
            #print coupled
            decoupled = tm.decoupled(coupled)
            print(self.str_matrix(decoupled))

    def tm_co_fitter(self, seeds):
        output = {}
        dim = len(seeds)
        tolerance = self.config_co["tolerance"]
        max_iter = self.config_co["max_iterations"]
        if max_iter == 0:
            return {
                "seed":seeds,
                "tm":[[0., 1., 0., 0., 0.],
                      [0., 0., 1., 0., 0.],
                      [0., 0., 0., 1., 0.],
                      [0., 0., 0., 0., 1.]],
                "substitutions":utils.utilities.do_lattice(self.config,
                                                           self.subs,
                                                           self.config_co["subs_overrides"]),
                "errors":[-1],
            }
        deltas = self.config_co["deltas"]
        a_track = None
        tm = None
        new_deltas = copy.deepcopy(deltas)
        errors = []
        new_seeds = seeds
        while len(errors) == 0 or (errors[-1] > tolerance and len(errors) < max_iter):
            print("----------------\nLooping with seed", self.str_matrix(new_seeds), end=' ')
            print("delta", self.str_matrix(deltas, "4.2g"), end=' ')
            print("error", self.get_error(new_deltas), "tolerance", tolerance)
            try:
                tm_list_in, tm_list_out, ref_track = self.get_tm(new_seeds, deltas)
            except RuntimeError:
                # if tm calculation fails, we abort with the last successful result
                # stored in seeds (e.g. tracking not convergent)
                break
            seeds = new_seeds 
            tm = self.fit_matrix(tm_list_in, tm_list_out)
            new_co, dm = self.get_co(tm)
            self.print_ref_track(ref_track, seeds, dm)
            for i in range(dim):
                new_deltas[i] = abs(tm_list_out[0][i] - tm_list_in[0][i])
                if self.config_co["adapt_deltas"] and new_deltas[i] < deltas[i]:
                    deltas[i] = new_deltas[i]
            errors.append(self.get_error(new_deltas))
            #new_co = self.check_stuck(errors, new_co, tm_list_in, tm_list_out, seeds)
            new_seeds = [seeds[i]+x for i, x in enumerate(new_co)]
        print("Finished iteration with deltas", deltas, "rms", sum([d*d for d in deltas])**0.5)
        if tm:
            tm = tm.get_coefficients_as_matrix()
        output = {
            "seed":seeds,
            "tm":tm,
            "substitutions":self.subs_tracking,
            "errors":errors,
        }
        return output

    def get_new_seed(self, config_seed):
        if len(self.output_list) == 0:
            return config_seed
        elif len(self.output_list) == 1:
            seed = Hit.new_from_dict(self.output_list[-1]["seed_hit"])
            return [seed[var] for var in self.var_list]
        else:
            seed0 = Hit.new_from_dict(self.output_list[-2]["seed_hit"])
            seed1 = Hit.new_from_dict(self.output_list[-1]["seed_hit"])
            s0 = [seed0[var] for var in self.var_list]
            s1 = [seed1[var] for var in self.var_list]
            s2 = [2*s1[i]-s0[i] for i in range(4)]
            return s2

    def save_track_orbit(self):
        file_name = self.config_co["orbit_file"]
        if file_name == "" or type(file_name) != type(""):
            return
        first, last = file_name.split(".")
        file_out = first+"_"+str(self.run_index)+"."+last
        os.rename(self.run_dir+"/"+file_name, self.run_dir+"/"+file_out)
        self.run_index += 1

    def get_subs_string(self):
        subs_dict = [{"subs":sub} for sub in self.config.substitution_list]
        subs_axes = utils.utilities.get_substitutions_axis(subs_dict, "subs")
        self.print_list = []
        for i, subs in enumerate(self.config.substitution_list):
            print_str = ""
            for axis in subs_axes:
                print_str += "    "+utils.utilities.sub_to_name(axis)
                print_str += utils.utilities.sub_to_units(axis)+": "
                print_str += str(subs_axes[axis][i])+"\n"
            self.print_list.append(print_str)

    def find_closed_orbits(self):
        output_uber_list = []
        self.get_subs_string()
        for config_seed in self.config_co["seed"]:
            self.output_list = []
            for i, self.subs in enumerate(self.config.substitution_list):
                print("\n\nNew closed orbit loop", i+1, "/", len(self.config.substitution_list), "with lattice values")
                print(self.print_list[i])
                self.energy = self.subs["__energy__"]
                self.get_tracking(False)
                seed = self.get_new_seed(config_seed)
                a_track = None
                try:
                    output = self.tm_co_fitter(seed)
                    output["seed_in"] = seed
                    output["seed_hit"] = self.seed_to_hit(output["seed"], 0.).dict_from_hit()
                    self.output_list.append(output)
                    a_track = self.track_one(output["seed"], 0., True)
                    output["ref_track"] = [hit.dict_from_hit() for hit in a_track]
                    self.print_ref_track(a_track, output["seed"], None)
                except Exception:
                    sys.excepthook(*sys.exc_info())
                self.save_track_orbit()
                self.save_output(self.output_list, False)
            output_uber_list += self.output_list
        self.save_output(output_uber_list, self.config_co["overwrite"])

    def save_output(self, output_list, overwrite):
        print("Overwriting closed orbits")
        tmp = self.run_dir+"/"+self.config_co["output_file"]+"_tmp"
        fout = open(tmp, "w")
        for output in output_list:
            print(json.dumps(output), file=fout)
        fout.close()
        if overwrite:
            output = self.config.run_control["output_dir"]+"/"+self.config_co["output_file"]
            os.rename(tmp, output)

    run_index = 1

def main(config):
    co_finder = ClosedOrbitFinder4D(config)
    co_finder.find_closed_orbits()

if __name__ == "__main__":
    main()
    if len(sys.argv) == 1:
        input()

