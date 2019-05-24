import glob
import sys
import os
import json
import time

import xboa.common as common
from xboa.hit import Hit
from opal_tracking import OpalTracking
from utils import utilities

class DAFinder(object):
    def __init__(self, config):
        self.closed_orbit_file_name = os.path.join(config.run_control["output_dir"], config.find_closed_orbits["output_file"])+".out"
        self.da_file_name = os.path.join(config.run_control["output_dir"], config.find_da["get_output_file"])
        self.scan_file_name = os.path.join(config.run_control["output_dir"], config.find_da["scan_output_file"]) 
        self.config = config
        self.run_dir =  os.path.join(config.run_control["output_dir"], config.find_da["run_dir"])
        self.co_list = self.load_closed_orbits()
        self.ref_hit = None
        self.min_delta = config.find_da["min_delta"]
        self.max_delta = config.find_da["max_delta"]
        self.data = []
        self.fout_scan_tmp = None
        self.fout_get_tmp = None
        self.required_n_hits = config.find_da["required_n_hits"]
        self.max_iterations = config.find_da["max_iterations"]
        self.tracking = self.setup()

    def get_all_da(self, co_index_list, seed_x, seed_y):
        if co_index_list == None:
            co_index_list = range(len(self.co_list))
        for i in co_index_list:
            try:
                co_element = self.co_list[i]
                print "Finding da for element", i
            except KeyError:
                print "Failed to find index", i, "in co_list of length", len(co_list)
                continue
            if seed_x != None and seed_x > 0.:
                co_element['x_da'] = self.get_da(co_element, 'x', seed_x)
            if seed_y != None and seed_y > 0.:
                co_element['y_da'] = self.get_da(co_element, 'y', seed_y)
            print >> self.fout_get(), json.dumps(co_element)
            self.fout_get().flush()

    def da_all_scan(self, co_index_list, x_list, y_list):
        if co_index_list == None:
            co_index_list = range(len(self.co_list))
        for i in co_index_list:
            try:
                co_element = self.co_list[i]
                print "Scanning da for element", i
            except KeyError:
                print "Failed to find index", i, "in co_list of length", len(co_list)
                continue
            co_element['da_scan'] = self.da_scan(co_element, x_list, y_list)

    def load_closed_orbits(self):
        fin = open(self.closed_orbit_file_name)
        co_list = [json.loads(line) for line in fin.readlines()]
        print "Loaded", len(co_list), "closed orbits"
        return co_list

    def setup(self):
        self.tmp_dir = "./"
        try:
            os.makedirs(self.run_dir)
        except OSError: # maybe the dir already exists
            pass
        os.chdir(self.run_dir)
        print "Running in", os.getcwd()
        self.opal_exe = os.path.expandvars("${OPAL_EXE_PATH}/opal")

    def reference(self, hit_dict):
        """
        Generate a reference particle
        """
        hit = Hit.new_from_dict(hit_dict)
        hit["x"] = 0.
        hit["px"] = 0.
        return hit

    def setup_tracking(self, co_element):
        subs = co_element["substitutions"]
        for item, key in self.config.find_da["subs_overrides"].iteritems():
            subs[item] = key
        print "Set up tracking for da with", 
        for key in sorted(subs.keys()):
            print utilities.sub_to_name(key), subs[key],
        print
        self.ref_hit = self.reference(co_element["hits"][0])
        lattice_src = self.config.tracking["lattice_file"]
        common.substitute(
            lattice_src, 
            self.run_dir+"/SectorFFAGMagnet.tmp",
            subs
        )
        tracking_file = self.config.find_da["probe_files"]
        self.tracking = OpalTracking(self.run_dir+"/SectorFFAGMagnet.tmp", self.tmp_dir+'/disttest.dat', self.ref_hit, tracking_file, self.opal_exe, self.tmp_dir+"/log")

    def new_seed(self):
        if self.data[-1][0] < self.min_delta: # reference run?
            return None
        if self.test_pass(*self.data[-1]): # upper limit is okay; keep going up
            if self.data[-1][0] > self.max_delta: # too big; give up
                return None
            return self.data[-1][0]*2.
        elif not self.test_pass(*self.data[0]): # lower limit is bad; try going down
            if abs(self.data[0][0]) < self.min_delta:
                return None
            return self.data[0][0]/2.
        else:
            for i, item in enumerate(self.data[1:]):
                if not self.test_pass(*item):
                    break
            if abs(self.data[i][0]-self.data[i+1][0]) < self.min_delta:
                return None
            return (self.data[i][0]+self.data[i+1][0])/2.

    def test_pass(self, seed, hits_list):
        return len(hits_list) > self.required_n_hits

    def events_generator(self, co_element, x_list, y_list):
        co_hit = co_element["hits"][0]
        for x in x_list:
            for y in y_list:
                a_hit = self.ref_hit.deepcopy()
                a_hit['x'] = co_hit['x'] + x
                a_hit['px'] = co_hit['px']
                a_hit['y'] += y
                yield {"x":x, "y":y}, a_hit

    def da_scan(self, co_element, x_list, y_list):
        self.setup_tracking(co_element)
        gen = self.events_generator(co_element, x_list, y_list)
        self.data = []
        finished = False
        while not finished:
            event_list = []
            track_list = []
            try:
                while len(event_list) < 1:
                    track, event = gen.next()
                    track_list.append(track)
                    event_list.append(event)
            except StopIteration:
                finished = True
            if len(event_list) == 0:
                break
            many_tracks = self.tracking.track_many(event_list)
            for i, hits in enumerate(many_tracks):
                print "Tracked", len(hits), "total hits with track", track_list[i], "first event x px y py", hits[0]["x"], hits[0]["px"], hits[0]["y"], hits[0]["py"]
                self.data.append([track_list[i], [a_hit.dict_from_hit() for a_hit in hits]])
        print >> self.fout_scan(), json.dumps(self.data)
        self.fout_scan().flush()

    def get_da(self, co_element, axis, seed_x):
        is_ref = abs(seed_x) < 1e-6
        self.setup_tracking(co_element)
        self.data = []
        co_delta = {"x":0, "y":0}
        iteration = 0
        while seed_x != None and iteration < self.max_iterations:
            co_delta[axis] = seed_x
            my_time = time.time()
            a_hit = Hit.new_from_dict(co_element["hits"][0])
            a_hit[axis] += seed_x
            try:
                hits = self.tracking.track_one(a_hit)
            except RuntimeError:
                sys.excepthook(*sys.exc_info())
                print "Never mind, keep on going..."
            self.data.append([co_delta[axis], [a_hit.dict_from_hit() for a_hit in hits]])
            self.data = sorted(self.data)
            print "Axis", axis, "Seed", seed_x, "Number of cells hit", len(hits), "in", time.time() - my_time, "[s]"
            sys.stdout.flush()
            seed_x = self.new_seed()
            if is_ref:
                seed_x = None
            iteration += 1
        self.data = [list(item) for item in self.data]
        return self.data

    def fout_scan(self):
        if self.fout_scan_tmp == None:
            file_name = self.scan_file_name+".tmp"
            self.fout_scan_tmp = open(file_name, "w")
            print "Opened file", file_name
        return self.fout_scan_tmp

    def fout_get(self):
        if self.fout_get_tmp == None:
            file_name = self.da_file_name+".tmp"
            self.fout_get_tmp = open(file_name, "w")
            print "Opened file", file_name
        return self.fout_get_tmp

def main(config):
    finder = DAFinder(config)
    finder.da_all_scan(config.find_da["row_list"], config.find_da["scan_x_list"], config.find_da["scan_y_list"])
    finder.get_all_da(config.find_da["row_list"], config.find_da["x_seed"], config.find_da["y_seed"], )
    print "Done find da"
    sys.stdout.flush()
    return
