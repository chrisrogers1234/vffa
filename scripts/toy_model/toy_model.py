import argparse
import glob
import json
import bisect
import copy
import os
import time
import math
import subprocess

import scipy
import matplotlib
import numpy
import xboa.bunch
import xboa.common

from utils import decoupled_transfer_matrix 
decoupled_transfer_matrix.DecoupledTransferMatrix.det_tolerance=1
from utils import utilities
from foil_model.material import Material
from foil_model.particle import Particle

numpy.set_printoptions(edgeitems=30, linewidth=200)
P_MASS = 938.27208816


class ToyModel(object):
    def __init__(self):
        n_pulses = 1

        #Run 101 with args [0.2116535191136668, -1.1802958558145478, 1.5605421957991457]Y MAX 10.0
        # injection setup
        self.momentum = 75.
        self.max_turn = 50 #n_pulses*1.1
        self.number_pulses = n_pulses
        self.number_per_pulse = 1
        self.turn_bumper_index = [range(n_pulses), range(n_pulses)] # first list is the turn number, second list is the index of the bumper magnet setting
        angle_u = -0.11628784278921322
        angle_v = 1.4012763546351992
        if n_pulses == 1:
            self.default_bumps = [[angle_u, 1.3, angle_v, 1.3]] # position of the proton orbit in aa coordinates
        else:
            self.default_bumps = [[angle_u, (n_pulses-1.-i)*1.3/(n_pulses-1.), angle_v, i*1.3/(n_pulses-1)] for i in range(n_pulses)] # position of the proton orbit in aa coordinates
        for line in self.default_bumps:
            print(line)
        self.default_injection = [math.pi/2, 0., 0.0, 0.0] # position of the injection orbit in aa coordinates
        self.foil_edge = -1
        self.foil_angle = math.degrees(1.1362432178340065)

        self.m_index = 1.4 # m^-1
        self.pulse_emittance = 1e-9 #0.026
        self.dp_over_p = 0.0 #0.0013 # 0.00 # sigma
        self.max_dp_over_p = self.dp_over_p*3
        self.foil_material = "carbon"
        self.foil_column_density = 0.0 #5e-6 # g/cm^3
        self.pid = 2212
        self.n_foil_sigma = 3 # number of sigma; set to < 0 to disable foil adjustment

        # execution parameters
        self.verbose = 8
        self.n_cells = 10
        self.accumulate = True
        self.do_plots = True
        self.plot_frequency = 1 #self.max_turn # controls whether to plot frames
        self.do_movie = False
        self.sleep_time = 0.1
        self.f_size = 20
        self.l_size = 14
        self.momentum_offset_injection = False
        self.do_scattering = True
        self.do_energy_loss = True

        # internal data
        self.bump_fields = []
        self.bump_orbits = []
        self.bump_tms = []
        self.beam_data = None
        self.injection_orbit = None
        self.dp_over_p_data = []
        self.beam_injection_turn = []
        self.turn = 0
        self.foil_hits = [] # number of foil hits per particle
        self.foil_hit_positions = [] # positions of foil hits
        self.first_turn_positions = [] # positions of foil hits
        self.output_dir = ""
        self.setup_material()
        self.setup_subplots()

    def do_config(self, config):
        for key in config:
            if key in self.__dict__:
                self.__dict__[key] = config[key]
            else:
                raise(KeyError("Did not recognise config key '"+str(key)+"'"))
        self.setup_material()

    def foil_edge_from_beam(self, n_sigma):
        ellipse = self.injection_ellipse()
        ellipse = numpy.array([
                [ellipse[0][0], ellipse[2][0],],
                [ellipse[2][0], ellipse[2][2],],
            ])
        theta = math.radians(self.foil_angle)
        R = numpy.array([
                [math.cos(theta), math.sin(theta)],
                [-math.sin(theta), math.cos(theta)]
            ])
        ellipse = numpy.dot(R, ellipse)
        ellipse = numpy.dot(ellipse, R.transpose())
        self.foil_edge = max(ellipse[1][1]**0.5*n_sigma, 0.1)

    def setup_material(self):
        self.foil = Material()
        self.foil.set_material(self.foil_material)
        self.foil_thickness = self.foil_column_density/self.foil.density

    def setup_subplots(self):
        self.fig = matplotlib.pyplot.figure(figsize=(20, 10))
        self.fig2 = matplotlib.pyplot.figure(figsize=(20, 10))
        self.axes2 = [
            self.fig2.add_subplot(2, 3, 1),
            self.fig2.add_subplot(2, 3, 2),
            self.fig2.add_subplot(2, 3, 3),
            self.fig2.add_subplot(2, 3, 4),
            self.fig2.add_subplot(2, 3, 5),
            self.fig2.add_subplot(2, 3, 6),
        ]
        self.axes = [
            self.fig.add_subplot(2, 3, 1),
            self.fig.add_subplot(2, 6, 7),
            self.fig.add_subplot(2, 6, 8),
            self.fig.add_subplot(2, 3, 2),
            self.fig.add_subplot(2, 6, 9),
            self.fig.add_subplot(2, 6, 10),
            self.fig.add_subplot(2, 3, 3),
            self.fig.add_subplot(2, 6, 11),
            self.fig.add_subplot(2, 6, 12),
            None,
            None,
        ]
        self.axes[9] = self.axes[1].twinx()
        self.axes[10] = self.axes[2].twinx()
        self.axes[0].set_title('Real Space')
        self.axes[3].set_title('Decoupled Space')
        self.axes[6].set_title('Action-Angle Coordinates')

    def parse_args(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("--output_dir", type=str, default="")
        parser.add_argument("--bump_settings", type=str, default="")
        parser.add_argument("--closed_orbit", type=str)
        parser.add_argument("--foil_station", type=int, default=0)
        args = parser.parse_args()
        self.settings = args.bump_settings
        self.closed_orbit = args.closed_orbit
        self.foil_station = args.foil_station
        self.output_dir = args.output_dir+"/"

    def initialise(self):
        self.load_closed_orbit()
        self.load_bump_settings()
        self.foil_edge_from_beam(self.n_foil_sigma)
        if self.verbose > 4:
            tm = self.bump_tms[0]
            print("Using transfer map with determinant", numpy.linalg.det(tm.m), " m:")
            print(numpy.real(tm.m))
            print("t:")
            print(numpy.real(tm.t))
            print(numpy.linalg.det(tm.t[0:2, 0:2]), numpy.linalg.det(tm.t[2:4, 2:4]), numpy.linalg.det(tm.t),)
            print("r:")
            print(numpy.real(tm.r))
            print(numpy.linalg.det(tm.r))
            print("v_t:")
            print(numpy.real(tm.v_t))
            print(numpy.linalg.det(tm.v_t[0:2, 0:2]), numpy.linalg.det(tm.v_t[2:4, 2:4]), numpy.linalg.det(tm.v_t))
            print("Closed orbits:")
            for i, turn in enumerate(self.turn_bumper_index[1]):
                orbit = self.bump_orbits[turn]
                o_str = ' '.join([format(x, "8.4g") for x in orbit])
                aa_str = ' '.join([format(x, "8.4g") for x in tm.coupled_to_action_angle(orbit-self.injection_orbit)])
                dec_str = ' '.join([format(x, "8.4g") for x in tm.decoupled(orbit-self.injection_orbit)])
                print("  Orbit", i, o_str, "  action angle", aa_str, "  decoupled", dec_str)
        self.beam_data = None
        self.turn = 0
        self.beam_injection_turn = []
        self.foil_hits = []
        self.foil_hit_positions = []
        self.first_turn_positions = []
        self.dp_over_p_data = []

    def load_closed_orbit(self):
        if self.verbose > 3:
            print("Loading closed orbits", self.closed_orbit)
        closed_orbit = json.loads(open(self.closed_orbit).readline())
        tm = closed_orbit[0]["tm"]
        one_cell = numpy.array([row[1:5] for row in tm])
        tm = copy.deepcopy(one_cell)
        for i in range(self.n_cells-1):
            tm = numpy.dot(tm, one_cell)
        tm = decoupled_transfer_matrix.DecoupledTransferMatrix(tm)
        tm.r/=numpy.linalg.det(tm.r)**0.25
        tm.r_inv/=numpy.linalg.det(tm.r_inv)**0.25
        self.bump_tms  = [tm]*(max(self.turn_bumper_index[1])+1)

    def default_bump_settings(self):
        if self.verbose > 3:
            print("Using default bump settings")
        self.bump_orbits  = [None]*(max(self.turn_bumper_index[1])+1)
        for i, bump_index in enumerate(self.turn_bumper_index[1]):
            tm = self.bump_tms[bump_index]
            action_angle_bump = self.default_bumps[i]
            action_angle_bump[1] /= self.momentum
            action_angle_bump[3] /= self.momentum
            print("Default bump settings", action_angle_bump)
            coupled_bump = tm.action_angle_to_coupled(action_angle_bump)
            print("Default bump settings invert", tm.coupled_to_action_angle(coupled_bump))
            self.bump_orbits[bump_index] = coupled_bump
        self.injection_orbit = tm.action_angle_to_coupled(self.default_injection)

    def load_bump_settings(self):
        if self.settings == "":
            self.default_bump_settings()
            return
        file_list = sorted(glob.glob(self.settings))
        if self.verbose > 3:
            print("Loading", len(file_list), "bump settings",
                  file_list[0], "...", file_list[-1], end=" ... ")
        for my_file in file_list:
            line_in = open(my_file).readlines()[-1]
            line_json = json.loads(line_in)
            if line_json["optimisation_stage"]["optimisation_stage"] != 1:
                continue
            for hit in line_json["tracking"]:
                if hit[0] == self.foil_station:
                    self.bump_fields.append(line_json["bump_fields"])
                    self.bump_orbits.append(numpy.array(hit[1:]))
                    break
        self.injection_orbit = self.bump_orbits[0]
        if self.verbose > 3:
            print("Done")

    def run_toy_model(self):
        self.setup_output_dir()
        try:
            if self.verbose > 3:
                print("Running")
            while True:
                self.do_one_turn()
        except StopIteration:
            pass

    def do_one_turn(self):
        self.inject_one() # inject on the injection orbit
        if self.verbose > 7:
            print("Turn", self.turn)
            print("  Bump orbit", self.get_bump_orbit())
            print("  injection ", self.injection_orbit)
        if self.verbose > 10:
            print(self.beam_data)
        #self.check_foil() # check for foil hits and apply foil physics
        self.plot_one_turn() # plot data
        self.turn_one()  # do one turn transfer matrix phase advance
        self.update_turn() # add one turn to turn (index)

    def print_aa(self):
        aa_hits = self.convert_to_aa(self.beam_data)
        for hit in aa_hits:
            print("Print aa", hit)


    def setup_output_dir(self):
        if self.output_dir == "":
            self.output_dir = os.path.split(self.closed_orbit)[0]+"/toy_model/"
        utilities.clear_dir(self.output_dir)

    def injection_ellipse(self):
        return self.get_transfer_matrix().get_v_m([self.pulse_emittance/self.momentum]*2)

    def inject_one(self):
        number = self.number_per_pulse
        if self.pulse_emittance > 0.:
            events = numpy.zeros((number, 4))
            dp_over_p_list = [0 for i in range(number)]
            n_successes = 0
            ellipse = self.injection_ellipse()
            ellipse_inv = numpy.linalg.inv(ellipse)
            while n_successes < self.number_per_pulse:
                if self.dp_over_p > 1e-9:
                    dp_over_p = numpy.random.normal(0., self.dp_over_p)
                    if dp_over_p > self.max_dp_over_p:
                        continue
                else:
                    dp_over_p = 0.
                mean = copy.deepcopy(self.injection_orbit)
                if self.momentum_offset_injection:
                    mean[2] += dp_over_p/self.m_index*1000. # mm
                dp_over_p_list[n_successes] = dp_over_p
                events[n_successes] = numpy.random.multivariate_normal(mean, ellipse, 1)
                n_successes += 1
        else:
            raise ValueError("toy_model.pulse_emittance must be > 0")
            inject = copy.deepcopy(self.injection_orbit).tolist()
            events = numpy.array([inject for i in range(number)])
            dp_over_p_list = [0 for i in range(number)]
        if type(self.beam_data) == type(None):
            inject = copy.deepcopy(self.injection_orbit).tolist()
            self.beam_data = events
            self.dp_over_p_data = dp_over_p_list
            self.beam_injection_turn += [self.turn for i in range(number)]
            return
        if len(self.beam_data) >= self.number_pulses*self.number_per_pulse:
            return
        self.beam_data = numpy.append(self.beam_data,
                                      events,
                                      axis=0)
        self.dp_over_p_data += dp_over_p_list
        self.beam_injection_turn += [self.turn for i in range(number)]
        self.first_turn_positions += [(hit[0], hit[2]) for hit in events]

    def transform_to_relative_coordinates(self, beam_data = None):
        # transform to coordinate system of bumped orbit
        if type(beam_data) == type(None):
            beam_data = self.beam_data
        for i, row in enumerate(beam_data):
            bump_orbit = self.get_bump_orbit(self.turn, self.dp_over_p_data[i])
            row -= bump_orbit
        return beam_data

    def transform_to_absolute_coordinates(self, beam_data = None):
        # transform back to global coordinate system
        if type(beam_data) == type(None):
            beam_data = self.beam_data
        for i, row in enumerate(beam_data):
            bump_orbit = self.get_bump_orbit(self.turn, self.dp_over_p_data[i])
            row += bump_orbit
        return beam_data

    def turn_one(self):
        tm = numpy.real(self.get_transfer_matrix().m)
        #self.transform_to_relative_coordinates() # transform to coordinates relative to momentum-dependent closed orbit
        self.beam_data = numpy.transpose(self.beam_data)
        self.beam_data = numpy.dot(tm, self.beam_data)
        self.beam_data = numpy.transpose(self.beam_data)
        #self.transform_to_absolute_coordinates() # transform to global coordinate system

    def update_turn(self):
        if self.last_turn():
            raise StopIteration("turn "+str(self.turn)+" > max turn "+str(self.max_turn))
        self.turn += 1

    def last_turn(self):
        return self.turn >= self.max_turn

    def get_bump_orbit(self, turn = None, dp_over_p = 0):
        index = self.get_bumper_index(turn)
        bumper_1 = self.turn_bumper_index[1][int(index)]
        if int(index)+1 ==  len(self.turn_bumper_index[1]):
            return self.bump_orbits[bumper_1]
        bumper_2 = self.turn_bumper_index[1][int(index)+1]
        orbit_1 = self.bump_orbits[bumper_1]
        orbit_2 = self.bump_orbits[bumper_2]
        orbit = (orbit_2-orbit_1)*(index-int(index))+orbit_1
        orbit[2] += dp_over_p/self.m_index*1000.
        return orbit

    def get_bumper_index(self, turn = None):
        if turn == None:
            turn = self.turn
        index = bisect.bisect_left(self.turn_bumper_index[0], turn)
        if index <= 0:
            return 0
        elif index >= len(self.turn_bumper_index[0]):
            return len(self.turn_bumper_index[0])-1
        turn_0 = self.turn_bumper_index[0][index-1]
        turn_1 = self.turn_bumper_index[0][index]
        index = (turn-turn_0)/float(turn_1-turn_0)+self.turn_bumper_index[1][index-1]
        return index

    def get_transfer_matrix(self, turn = None):
        index = int(self.get_bumper_index(turn))
        bumper_setting = self.turn_bumper_index[1][index]
        return self.bump_tms[bumper_setting]

    def make_2d_hist(self, ax1, ax2, lab1, lab2, axes, hits, centre, range_):
        x_list = [hit[ax1] for hit in hits]
        y_list = [hit[ax2] for hit in hits]
        axes.clear()
        n_bins = int(max(len(x_list)**0.5/3, 10))
        axes.hist2d(x_list, y_list, n_bins)
        color = matplotlib.cm.get_cmap()(0)
        axes.set_facecolor(color)
        axes.set_xlim(-range_[ax1]+centre[ax1], range_[ax1]+centre[ax1])
        axes.set_ylim(-range_[ax2]+centre[ax2], range_[ax2]+centre[ax2])
        axes.set_xlabel(lab1, fontsize=self.f_size)
        axes.set_ylabel(lab2, fontsize=self.f_size)
        axes.tick_params(labelsize = self.l_size)

    def make_a_plot(self, ax1, ax2, lab1, lab2, axes, hits, centre, range_):
        #add p1 off the edge of the plot; this forces the z axis
        p1 = [2*range_[ax1]+centre[ax1], 2*range_[ax2]+centre[ax2], 1.]
        const = [1, self.momentum, 1, self.momentum]
        x_list = [p1[0]]+[hit[ax1]*const[ax1] for hit in hits]
        y_list = [p1[1]]+[hit[ax2]*const[ax2] for hit in hits]
        z_list = [1.]+[t/self.number_pulses for t in self.beam_injection_turn]
        colors = matplotlib.pyplot.cm.coolwarm
        if not self.accumulate:
            axes.clear()
        marker_size = 100/len(hits)**0.5
        axes.scatter(x_list, y_list, c=z_list, marker="o", cmap=colors, s=marker_size)
        axes.set_xlim(-range_[ax1]+centre[ax1], range_[ax1]+centre[ax1])
        axes.set_ylim(-range_[ax2]+centre[ax2], range_[ax2]+centre[ax2])
        axes.set_xlabel(lab1, fontsize=self.f_size)
        axes.set_ylabel(lab2, fontsize=self.f_size)

    def trajectory_plot(self, ax1, ax2, axes, trajectory, specials):
        x_list = [hit[ax1] for hit in trajectory]
        y_list = [hit[ax2] for hit in trajectory]
        axes.plot(x_list, y_list)
        x_list = [hit[ax1] for hit in specials]
        y_list = [hit[ax2] for hit in specials]
        axes.plot(x_list, y_list, 'o', fillstyle='none')

    def make_1d_plot(self, ax1, lab1, axes, hits, centre, range_):
        x_list = [hit[ax1] for hit in hits]
        n_bins = int(max(10, len(x_list)/10))
        axes.clear()
        axes.frameon = False
        axes.get_xaxis().set_visible(False)
        axes.get_yaxis().set_visible(False)
        axes.hist(x_list, bins=n_bins)
        axes.set_xlim(-range_[ax1]+centre[ax1], range_[ax1]+centre[ax1])
        axes.set_xlabel(lab1, fontsize=self.f_size)
        y_max = len(x_list)/n_bins*10
        axes.set_ylim(0, y_max)
        print("Y MAX", y_max)

    def plot_bump_orbit(self):
        n = self.number_pulses
        turn = range(self.number_pulses)
        x_list, px_list, y_list, py_list = [None]*n, [None]*n, [None]*n, [None]*n
        u_list, pu_list, v_list, pv_list = [None]*n, [None]*n, [None]*n, [None]*n
        for i in turn:
            matrix = self.get_transfer_matrix(i)
            orbit = self.get_bump_orbit(i)
            x_list[i] = orbit[0]
            px_list[i] = orbit[1]
            y_list[i] = orbit[2]
            py_list[i] = orbit[3]
            decoupled_orbit = matrix.decoupled(orbit)
            u_list[i] = decoupled_orbit[0]
            pu_list[i] = decoupled_orbit[1]
            v_list[i] = decoupled_orbit[2]
            pv_list[i] = decoupled_orbit[3]
        figure = matplotlib.pyplot.figure(figsize=(20, 10))
        axes = figure.add_subplot(2, 2, 1)
        axes.plot(turn, x_list, label="x [mm]")
        axes.plot(turn, y_list, label="y [mm]")
        axes.set_xlabel("Turn", fontsize=self.f_size)
        axes.set_ylabel("", fontsize=self.f_size)
        axes.tick_params(labelsize = self.l_size)
        axes.legend()
        axes = figure.add_subplot(2, 2, 2)
        axes.plot(turn, px_list, label="$p_x$ [MeV/c]")
        axes.plot(turn, py_list, label="$p_y$ [MeV/c]")
        axes.set_xlabel("Turn", fontsize=self.f_size)
        axes.set_ylabel("", fontsize=self.f_size)
        axes.tick_params(labelsize = self.l_size)
        axes.legend()
        axes = figure.add_subplot(2, 2, 3)
        axes.plot(turn, u_list, label="u")
        axes.plot(turn, v_list, label="v")
        axes.set_xlabel("Turn", fontsize=self.f_size)
        axes.set_ylabel("", fontsize=self.f_size)
        axes.tick_params(labelsize = self.l_size)
        axes.legend()
        axes = figure.add_subplot(2, 2, 4)
        axes.plot(turn, pu_list, label="$p_u$")
        axes.plot(turn, pv_list, label="$p_v$")
        axes.set_xlabel("Turn", fontsize=self.f_size)
        axes.set_ylabel("", fontsize=self.f_size)
        axes.tick_params(labelsize = self.l_size)
        axes.legend()
        figure.savefig(self.output_dir+"turn_position.png")

    def plot_spherical(self, ax1, ax2, ax3, lab1, lab2, lab3):
        tm = self.get_transfer_matrix()
        print("Plot spherical")
        aa_hits = [tm.coupled_to_nd_action_angles(hit) for hit in self.beam_data]
        x_list = [math.degrees(hit[ax1]) for hit in aa_hits]
        y_list = [math.degrees(hit[ax2]) for hit in aa_hits]
        z_list = [math.degrees(hit[ax3]) for hit in aa_hits]
        colors = matplotlib.pyplot.cm.coolwarm
        axes = self.axes2[5]
        marker_size = 100/len(aa_hits)**0.5
        scat = axes.scatter(x_list, y_list, c=z_list, marker='o', cmap=colors, s=marker_size)
        axes.set_xlabel(lab1, fontsize=self.f_size)
        axes.set_ylabel(lab2, fontsize=self.f_size)
        axes.tick_params(labelsize = self.l_size)
        self.fig2.colorbar(scat, label=lab3)

        r_list = [hit[3]*self.momentum for hit in aa_hits]
        axes = self.axes2[2]
        n_bins = 50
        scat = axes.hist(r_list, n_bins)
        n_bins = max(int(n_bins*axes.get_ylim()[1]/100.), 10)
        axes.clear()
        scat = axes.hist(r_list, n_bins)
        axes.set_xlabel("$A_{4D}$ [$\\mu$m]", fontsize=self.f_size)
        axes.tick_params(labelsize = self.l_size)

    def amplitude_plots(self, aa_hits):
        #self.test_amplitude()
        #input("Press <CR> to continue")
        for hit in aa_hits:
            print("Plot amplitude", hit)
        aa_range = [200, 1.5, 200, 1.5]
        aa_centre = [0., aa_range[1]-0.1, 0., aa_range[3]-0.1]
        self.make_a_plot(1, 3, "$A_{u}$ [$\\mu$m]", "$A_{v}$ [$\\mu$m]", self.axes[6], aa_hits, aa_centre, aa_range)
        self.make_a_plot(0, 1, "$\\phi_{u}$ [deg]", "$A_{u}$ [$\\mu$m]", self.axes[7], aa_hits, aa_centre, aa_range)
        self.make_a_plot(2, 3, "$\\phi_{v}$ [deg]", "$A_{v}$ [$\\mu$m]", self.axes[8], aa_hits, aa_centre, aa_range)
        if self.verbose > 8:
            print("AA mean", numpy.mean(aa_hits, 0), "AA range", numpy.ptp(aa_hits, 0))

    def convert_to_decoupled(self, hits):
        r_mat = self.get_transfer_matrix().r_inv
        hits = numpy.dot(r_mat, numpy.transpose(hits))
        hits = numpy.transpose(numpy.real(hits))
        return hits

    def convert_to_aa(self, coupled_hits):
        norm = [180./math.pi, P_MASS/1000., 180./math.pi, P_MASS/1000.] # normalised action in microns
        tm = self.get_transfer_matrix()
        aa_hits = [tm.coupled_to_action_angle(hit) for hit in coupled_hits]
        aa_hits = [[u*norm[i] for i, u in enumerate(hit)] for hit in aa_hits]
        return aa_hits

    def draw_foil(self, axes):
        dw = 100.0
        dt = self.foil_edge
        theta = math.radians(self.foil_angle)
        x = [dt*math.sin(theta)+dw*math.cos(theta), dt*math.sin(theta)-dw*math.cos(theta)]
        y = [dt*math.cos(theta)-dw*math.sin(theta), dt*math.cos(theta)+dw*math.sin(theta)]
        axes.plot(x, y, color="gray")

    def plot_one_turn(self):
        if not self.do_plots:
            return
        if not self.last_turn() and self.turn % self.plot_frequency != 0:
            return

        hits = copy.deepcopy(self.beam_data)
        real_range = self.real_range
        real_centre = self.real_centre

        dec_range = self.dec_range
        dec_centre = self.dec_centre

        self.fig.suptitle("Turn "+str(self.turn)+"    Bumper setting "+str(self.get_bumper_index()+1))
        if self.verbose > 8:
            print("Real mean", numpy.mean(hits, 0), "Real range", numpy.ptp(hits, 0))
        trajectory = numpy.array([self.bump_orbits[i] for i in self.turn_bumper_index[1]])
        specials = numpy.array([self.injection_orbit, self.get_bump_orbit()])
        self.make_2d_hist(0, 2, "x [mm]", "y [mm]", self.axes2[0], hits, real_centre, real_range)
        self.make_a_plot(0, 2, "x [mm]", "y [mm]", self.axes[0], hits, real_centre, real_range)
        self.trajectory_plot(0, 2, self.axes[0], trajectory, specials)
        self.draw_foil(self.axes[0])
        self.make_a_plot(0, 1, "x [mm]", "Conjugate momentum [MeV/c]", self.axes[1], hits, real_centre, real_range)
        self.make_1d_plot(0, "x [mm]", self.axes[9], hits, real_centre, real_range)
        self.trajectory_plot(0, 1, self.axes[1], trajectory, specials)
        self.make_a_plot(2, 3, "y [mm]", "", self.axes[2], hits, real_centre, real_range)
        self.make_1d_plot(2, "y [mm]", self.axes[10], hits, real_centre, real_range)
        self.trajectory_plot(2, 3, self.axes[2], trajectory, specials)

        relative_hits = self.transform_to_relative_coordinates(copy.deepcopy(self.beam_data))

        decoupled_hits = self.convert_to_decoupled(relative_hits)
        trajectory = self.convert_to_decoupled(trajectory) # should be relative coordinates?
        specials = self.convert_to_decoupled(specials) # should be relative coordinates?

        if self.verbose > 8:
            print("Decoupled mean", numpy.mean(decoupled_hits, 0), "Decoupled range", numpy.ptp(decoupled_hits, 0))
        self.make_a_plot(0, 2, "u", "v", self.axes[3], decoupled_hits, dec_centre, dec_range)
        self.trajectory_plot(0, 2, self.axes[3], trajectory, specials)
        tune_u = self.get_transfer_matrix().get_phase_advance(0)/math.pi/2.
        tune_v = self.get_transfer_matrix().get_phase_advance(1)/math.pi/2.
        self.make_a_plot(0, 1, "u", "", self.axes[4], decoupled_hits, dec_centre, dec_range)
        self.trajectory_plot(0, 1, self.axes[4], trajectory, specials)
        self.axes[4].text(0.05, 0.95, "$\\nu_u$: "+format(tune_u, "6.4g"), transform=self.axes[4].transAxes)
        self.make_a_plot(2, 3, "v", "", self.axes[5], decoupled_hits, dec_centre, dec_range)
        self.trajectory_plot(2, 3, self.axes[5], trajectory, specials)
        self.axes[5].text(0.05, 0.95, "$\\nu_v$: "+format(tune_v, "6.4g"), transform=self.axes[5].transAxes)
        self.print_aa()

        aa_hits = self.convert_to_aa(self.beam_data)
        self.amplitude_plots(aa_hits)

        self.transform_to_absolute_coordinates()

        name = "turn_"+str(self.turn).rjust(5, "0")+".png"
        self.fig.savefig(self.output_dir+name)
        if self.sleep_time > 0:
            matplotlib.pyplot.pause(self.sleep_time)

    def plot_foil_hits(self, range_, centre):
        n_bins = 50 #max(self.foil_hits)+2
        axes = self.axes2[1]
        axes.hist(self.foil_hits, n_bins, (-0.5, n_bins-0.5), rwidth=0.2)
        axes.set_xlabel("Number of foil hits", fontsize=self.f_size)
        axes.tick_params(labelsize = self.l_size)
        axes = self.axes2[4]
        axes.hist(self.dp_over_p_data, n_bins, (-0.02, 0.02))
        axes.set_xlabel("dp/p", fontsize=self.f_size)
        axes.tick_params(labelsize = self.l_size)
        axes = self.axes2[3]
        x_list = [pos[0] for pos in self.foil_hit_positions]
        y_list = [pos[1] for pos in self.foil_hit_positions]
        marker_size = 10
        if len(self.foil_hit_positions) > 0:
            marker_size = 100/len(self.foil_hit_positions)**0.5
        axes.scatter(x_list, y_list, marker="o", s=marker_size)

        x_list = [pos[0] for pos in self.first_turn_positions]
        y_list = [pos[1] for pos in self.first_turn_positions]
        axes.scatter(x_list, y_list, marker="o", s=marker_size)

        self.draw_foil(axes)

        axes.set_xlim(-range_[0]+centre[0], range_[0]+centre[0])
        axes.set_ylim(-range_[2]+centre[2], range_[2]+centre[2])
        axes.set_xlabel("x [mm]", fontsize=self.f_size)
        axes.set_ylabel("y [mm]", fontsize=self.f_size)

    def foil_model(self, event_number):
        if not self.do_scattering and not self.do_energy_loss:
            return
        momentum = self.momentum*(1+self.dp_over_p_data[event_number])
        particle = Particle.new_from_momentum(momentum, self.pid)
        if self.do_scattering:
            sigma = self.foil.scattering(particle, self.foil_thickness)
            scat = numpy.random.randn(2)*sigma
            self.beam_data[event_number][1] += scat[0]
            self.beam_data[event_number][3] += scat[1]
            if self.verbose > 5 and event_number == 0:
                print("scattering sigma", sigma, "angles", scat, "psv", self.beam_data[event_number])
        if self.do_energy_loss:
            dE = self.foil.energy_loss_dz(particle)*self.foil_thickness #dE negative
            energy = particle.set_kinetic_energy(particle.get_kinetic_energy()+dE)
            dp_over_p = (particle.get_momentum()-self.momentum)/self.momentum
            if self.verbose > 5 and event_number == 0:
                print("Losing", dE, "of energy leaving", dp_over_p, "dp/p")
            self.dp_over_p_data[event_number] = dp_over_p

    def inside_foil(self, psv):
        theta = math.radians(self.foil_angle)
        dy = psv[2]*math.cos(theta)+psv[0]*math.sin(theta)
        return dy - self.foil_edge < 0

    def check_foil(self):
        self.foil_hits += [0]*(len(self.beam_data)-len(self.foil_hits))
        if self.verbose > 8:
            print("For foil", foil_dimensions)
        for i, hit in enumerate(self.beam_data):
            if self.verbose > 8:
                print("    Check foil", hit, end=" ")
            if self.inside_foil(hit):
                self.foil_hits[i] += 1
                self.foil_hit_positions.append((hit[0], hit[2]))
                self.foil_model(i)
                if self.verbose > 8:
                    print("hit")
            else:
                if self.verbose > 8:
                    print("miss")

    def movie(self):
        here = os.getcwd()
        os.chdir(self.output_dir)
        #mencoder mf://turn*.png -mf w=800:h=600:fps=5:type=png -ovc lavc -lavcopts vcodec=msmpeg4:mbd=2:trell -oac copy -o injection.avi
        try:
            output = subprocess.check_output(["mencoder",
                                    "mf://*.png",
                                    "-mf", "w=800:h=600:fps=5:type=png",
                                    "-ovc", "lavc",
                                    "-lavcopts", "vcodec=msmpeg4:vbitrate=2000:mbd=2:trell",
                                    "-oac", "copy",
                                    "-o", "injection.avi"])
        except:
            print("bob")
        os.chdir(here)

    def finalise(self):
        if not self.do_plots:
            return
        self.plot_foil_hits(self.real_range, self.real_centre)
        self.plot_spherical(0, 1, 2, "$\\phi_1$ [deg]", "$\\phi_2$  [deg]", "$\\phi_3$  [deg]")
        self.plot_bump_orbit()
        if self.max_turn/self.plot_frequency > 1 and self.do_movie:
            self.movie()
        self.fig.savefig(self.output_dir+"phase_space_final.png")
        self.fig2.savefig(self.output_dir+"summary_final.png")

    real_range = [10., 1, 20, 0.5]
    real_centre = [0.0, 0.0, 5.0, 0.0]

    dec_range = [5., 0.5, 15, 0.25]
    dec_centre = [0., 0., 0., 0.]



def main():
    model = ToyModel()
    model.parse_args()
    model.initialise()
    model.run_toy_model()
    model.finalise()
    matplotlib.pyplot.show(block=False)
    input("Press <CR> to finish")

if __name__ == "__main__":
    main()

