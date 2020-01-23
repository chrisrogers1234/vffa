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

#Jobs:
#* rotating planes might work; now need to add option to tracking to do analysis in reference coordinate system
#* fix RK4 - comparison is pretty good between RK4 and tracking, but not perfect. Check!

class ClosedOrbitFinder4D(object):
    def __init__(self, config):
        self.config = config
        self.config_co = config.find_closed_orbits
        self.energy = self.config.substitution_list[0]["__energy__"]
        self.run_dir = os.path.join(self.config.run_control["output_dir"],
                                    self.config_co["run_dir"])
        self.centroid = None
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

    @classmethod
    def rotation_matrix(cls, r1, r2):
        # rotation matrix that rotates r1 onto r2
        # note this is not uniquely defined
        v = numpy.cross(r1, r2)/numpy.linalg.norm(r1)/numpy.linalg.norm(r2)
        st = numpy.linalg.norm(v)
        v /= st
        ct = math.cos(math.asin(st))
        rotation_matrix = [
            [v[0]*v[0]*(1-ct)+ct,      -v[2]*st+v[0]*v[1]*(1-ct), v[0]*v[2]*(1-ct)+v[1]*st],
            [v[0]*v[1]*(1-ct)+v[2]*st,       ct+v[1]*v[1]*(1-ct), v[1]*v[2]*(1-ct)-v[0]*st],
            [v[0]*v[2]*(1-ct)-v[1]*st,  v[0]*st+v[1]*v[2]*(1-ct), v[2]*v[2]*(1-ct)+ct],
        ]
        return numpy.array(rotation_matrix)


    def rotate_from_centroid(self, hit_list):
        tm_list = []
        # 0 element in the hit_list defines the centroid
        # take s-vector as being the momentum vector of hit_list[0]
        # take all positions, etc perpendicular to the s-vector
        centroid = hit_list[0]
        ref_position = numpy.array([centroid['x'],
                                    centroid['y'],
                                    centroid['z']])
        ref_momentum = numpy.array([centroid['px'],
                                    centroid['py'],
                                    centroid['pz']])
        svec = ref_momentum/numpy.linalg.norm(ref_momentum)
        #rot = self.rotation_matrix(numpy.array([0., 0., 1.]), svec)
        rot = numpy.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
        for hit in hit_list:
            # project onto plane normal to the momentum vector
            hit_position = numpy.array([hit['x'], hit['y'], hit['z']])
            #distance_to_plane = numpy.dot(svec, hit_position-ref_position)
            #hit['x'] += distance_to_plane*hit['px']
            #hit['y'] += distance_to_plane*hit['py']
            #hit['z'] += distance_to_plane*hit['pz']
            # convert to coordinate system relative to centroid
            hit_pos = numpy.array([
                hit['x'],# - ref_position[0],
                hit['y'],# - ref_position[1],
                hit['z'],# - ref_position[2],
            ])
            hit_mom = numpy.array([
                hit['px'],
                hit['py'],
                hit['pz'],
            ])
            # rotate to coordinate system parallel to centroid
            hit_pos = numpy.dot(rot, hit_pos)
            hit_mom = numpy.dot(rot, hit_mom)
            vector = [hit_pos[0], hit_mom[0], hit_pos[1], hit_mom[1]]
            tm_list.append(vector)
        return tm_list, rot

    def track_many(self, seed_list, t, is_final):
        overrides = self.config_co["subs_overrides"]
        if is_final:
            overrides = self.config_co["final_subs_overrides"]
        overrides["__n_particles__"] = len(seed_list)+1
        hit_list = []
        for seed in seed_list:
            hit = self.seed_to_hit(seed, t)
            hit_list.append(hit)
        hit_list.insert(0, hit_list[0])
        os.chdir(self.run_dir)
        self.subs_tracking = utils.utilities.do_lattice(self.config, self.subs, overrides)
        track_list = self.tracking.track_many(hit_list)
        os.chdir(self.here)
        return track_list

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

    def get_values(self, dim): # makes a hyperdiamond at +- 1
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
            if do_break:
                break
            if do_continue:
                continue
            coefficients += [PolynomialCoefficient(sum_index, j, 0.0) for j in range(4)]
        tm = PolynomialMap.least_squares(tm_list_in, tm_list_out, 1, coefficients)
        return tm

    def get_tm(self, seeds, deltas):
        us_cell = self.config_co["us_cell"]
        ds_cell = self.config_co["ds_cell"]
        dim = len(seeds)
        tm_list_in, tm_list_out = [], []
        value_list = []
        for j, values in enumerate(self.get_values(dim)):
            value_list.append([seeds[i]+values[i]*deltas[i] for i in range(dim)])
        track_list = self.track_many(value_list, 0., False)[1:]
        try:
            tm_list_in, rot_in = self.rotate_from_centroid([a_track[us_cell] for a_track in track_list])
            tm_list_out, rot_out = self.rotate_from_centroid([a_track[ds_cell] for a_track in track_list])
            print(self.str_matrix(tm_list_in, "14.8g"))
            print(self.str_matrix(tm_list_out, "14.8g"))
        except IndexError:
            err = "Output had "+str(len(track_list))+" tracks with"
            err += str([len(track) for track in track_list])+" track points. "
            err += "Require "+str( ds_cell+1 )+" track points."
            print(err)
            sys.exit()
            raise RuntimeError(err) from None
        track_list[0]
        return tm_list_in, rot_in, tm_list_out, rot_out, track_list[0]

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

    def rotate_co(self, new_co, rot):
        pos = [new_co[0], new_co[2], 0.]
        return new_co

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
            self.centroid = self.seed_to_hit(new_seeds, 0.)
            try:
                print("tracking in")
                tm_list_in, rot_in, tm_list_out, rot_out, ref_track = self.get_tm(new_seeds, deltas)
                print("tracking out")
            except Exception:
                # if tm calculation fails, we abort with the last successful result
                # stored in seeds (e.g. tracking not convergent)
                sys.excepthook(*sys.exc_info())
                break
            seeds = new_seeds
            tm = self.fit_matrix(tm_list_in, tm_list_out)
            new_co, dm = self.get_co(tm)
            print("New closed orbit", new_co)
            #self.rotate_co(new_co, rot_in)
            self.print_ref_track(ref_track, seeds, dm)
            for i in range(dim):
                new_deltas[i] = abs(tm_list_out[0][i] - tm_list_in[0][i])
                if self.config_co["adapt_deltas"] and new_deltas[i] < deltas[i]:
                    deltas[i] = new_deltas[i]
            errors.append(self.get_error(new_deltas))
            #new_seeds = [seeds[i]+x for i, x in enumerate(new_co)]
            new_seeds = [x for i, x in enumerate(new_co)]
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
        if len(self.output_list) == 0 or len(self.output_list[-1]) == 0:
            return config_seed
        elif len(self.output_list[-1]) == 1:
            seed = Hit.new_from_dict(self.output_list[-1][-1]["seed_hit"])
            return [seed[var] for var in self.var_list]
        else:
            seed0 = Hit.new_from_dict(self.output_list[-1][-2]["seed_hit"])
            seed1 = Hit.new_from_dict(self.output_list[-1][-1]["seed_hit"])
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
        self.output_list = []
        self.get_subs_string()
        for config_seed in self.config_co["seed"]:
            self.output_list.append([])
            for i, self.subs in enumerate(self.config.substitution_list):
                print("\n\nNew closed orbit loop", i+1, "/", len(self.config.substitution_list), "with lattice values")
                print(self.print_list[i])
                self.energy = self.subs["__energy__"]
                self.get_tracking(False)
                seed = self.get_new_seed(config_seed)
                self.centroid = self.seed_to_hit(seed, 0.)
                a_track = None
                try:
                    output = self.tm_co_fitter(seed)
                    output["seed_in"] = seed
                    output["seed_hit"] = self.seed_to_hit(output["seed"], 0.).dict_from_hit()
                    self.output_list[-1].append(output)
                    a_track = self.track_many([output["seed"]]*3, 0., True)[1]
                    output["ref_track"] = [hit.dict_from_hit() for hit in a_track]
                    self.print_ref_track(a_track, output["seed"], None)
                except Exception:
                    sys.excepthook(*sys.exc_info())
                self.save_track_orbit()
                self.save_output(self.output_list, False)
        self.save_output(self.output_list, self.config_co["overwrite"])

    def save_output(self, output_list, overwrite):
        print("Overwriting closed orbits")
        output_dir = self.config.run_control["output_dir"]
        tmp = output_dir+"/"+self.config_co["output_file"]+"_tmp"
        fout = open(tmp, "w")
        for output in output_list:
            print(json.dumps(output), file=fout)
        print("Saved to", tmp)
        fout.close()
        if overwrite:
            output = output_dir+"/"+self.config_co["output_file"]
            os.rename(tmp, output)

    run_index = 1

def main(config):
    co_finder = ClosedOrbitFinder4D(config)
    co_finder.find_closed_orbits()

if __name__ == "__main__":
    main()
    if len(sys.argv) == 1:
        input()

