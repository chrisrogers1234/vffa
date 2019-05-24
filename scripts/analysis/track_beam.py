import os
import json
import subprocess
import numpy.random
import xboa.common
from xboa.hit import Hit

from opal_tracking import OpalTracking
from opal_tracking import TunesAnalysis
from opal_tracking import PhaseSpacePlots
from utils import utilities

class TrackBeam(object):
    def __init__(self, config):
        self.config = config
        self.output_dir = config.run_control["output_dir"]
        self.centre = None
        self.ellipse = None
        self.run_dir = ""
        self.cwd = os.getcwd()
        self.hits_in = []
        self.hits_out = []
        self.energy = None

    def load_tune_data(self):
        file_name = self.output_dir+"/"+self.config.find_tune["output_file"]
        fin = open(file_name)
        data = [json.loads(line) for line in fin.readlines()]
        return data

    def fit_tune_data(self, data):
        eps_max = self.config.track_beam['eps_max']
        x_emittance = self.config.track_beam['x_emittance']
        y_emittance = self.config.track_beam['y_emittance']
        sigma_pz = self.config.track_beam['sigma_pz']
        sigma_z = self.config.track_beam['sigma_z']

        self.energy = data['substitutions']['__energy__']
        pid = self.config.tracking["pdg_pid"]
        mass = xboa.common.pdg_pid_to_mass[abs(pid)]
        p = ((self.energy+mass)**2 - mass**2)**0.5

        x_centre, x_ellipse = xboa.common.fit_ellipse(
                                                  data['x_signal'],
                                                  eps_max,
                                                  verbose = False)
        y_centre, y_ellipse = xboa.common.fit_ellipse(
                                                  data['y_signal'],
                                                  eps_max,
                                                  verbose = False)
        x_ellipse *= (x_emittance/numpy.linalg.det(x_ellipse))**0.5
        y_ellipse *= (y_emittance/numpy.linalg.det(y_ellipse))**0.5
        self.centre = numpy.array([x for x in x_centre]+[y for y in y_centre]+[0., p])
        self.ellipse = numpy.zeros((6, 6))
        for i in range(2):
            for j in range(2):
                self.ellipse[i, j] = x_ellipse[i, j]
                self.ellipse[i+2, j+2] = y_ellipse[i, j]
        #for i in [0, 2]:
        #    self.ellipse[i+1, i+1] *= p*p
        #    self.ellipse[i, i+1] *= p
        #    self.ellipse[i+1, i] *= p
        self.ellipse[4, 4] = sigma_z
        self.ellipse[5, 5] = sigma_pz
        print "Centre", self.centre
        print "Ellipse"
        print self.ellipse

    def setup_workspace(self):
        self.run_dir = self.config.run_control["output_dir"]+"/"+self.config.track_beam["run_dir"]
        try:
            os.makedirs(self.run_dir)
        except OSError:
            pass # maybe the dir already exists
        os.chdir(self.run_dir)

    def generate_beam(self):
        n_events = self.config.track_beam["subs_overrides"]["__n_events__"]
        events = numpy.random.multivariate_normal(self.centre, self.ellipse, n_events)
        keys = "x", "px", "y", "py", "z", "pz"
        self.hits_in = []
        for item in events:
            hit = self.reference()
            for i, key in enumerate(keys):
                hit[key] = item[i]
            self.hits_in.append(hit)
        for hit in self.hits_in[0:10]:
            print "   ", [hit[key] for key in keys]
        print "Made", len(self.hits_in), "hits"
        

    def run_tracking(self, index):
        opal_exe = self.config.tracking["opal_path"]
        input_file = self.config.tracking["lattice_file"]
        n_cores = self.config.tracking["n_cores"]
        mpi_exe = self.config.tracking["mpi_exe"]
        lattice_file = self.run_dir+'SectorFFAGMagnet.tmp'

        subs = self.config.substitution_list[index]
        for key, value in self.config.track_beam["subs_overrides"].iteritems():
            subs[key] = value
        xboa.common.substitute(input_file, lattice_file, subs)
        log_name = self.run_dir+"/log"
        ref_hit = self.reference()
        probe_files = self.config.track_beam["probe_files"]
        self.tracking = OpalTracking("SectorFFAGMagnet.tmp", 'disttest.dat', ref_hit, probe_files, opal_exe, log_name, None, n_cores, mpi_exe)

        tunes_analysis = TunesAnalysis(self.config)
        phase_space_plots = PhaseSpacePlots(self.config)
        tunes_analysis.set_match(self.centre[0:4], self.ellipse[0:4, 0:4])
        if self.config.track_beam["do_track"]:
            print "Running tracking with\n   ",
            for key, value in subs.iteritems():
                print utilities.sub_to_name(key)+":", value,
            print
            self.tracking.track_many(self.hits_in, None)
        print os.getcwd(), probe_files
        self.tracking._read_probes(tunes_analysis)
        #self.tracking._read_probes(phase_space_plots)


    def reference(self):
        """
        Generate a reference particle
        """
        hit_dict = {}
        hit_dict["pid"] = self.config.tracking["pdg_pid"]
        hit_dict["mass"] = xboa.common.pdg_pid_to_mass[abs(hit_dict["pid"])]
        hit_dict["charge"] = 1
        hit_dict["x"] = 0.
        hit_dict["kinetic_energy"] = self.energy
        return Hit.new_from_dict(hit_dict, "pz")

    def track(self):
        try:
            data = self.load_tune_data()
            self.setup_workspace()
            for i, item in enumerate(data):
                self.fit_tune_data(item)
                self.generate_beam()
                self.run_tracking(i)
        except:
            raise
        finally:
            os.chdir(self.cwd)

def main(config):
    tracker = TrackBeam(config)
    tracker.track()
