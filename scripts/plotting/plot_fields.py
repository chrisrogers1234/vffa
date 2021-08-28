"""
Plot a single closed orbit (once tracking has finished)
"""
import sys
import copy
import os
import math
import argparse
import h5py
import glob
import shutil


import matplotlib
import ROOT

import xboa.common

import PyOpal.parser
import PyOpal.field
import utils.utilities

class GetFields(object):
    def __init__(self, lattice_file):
        self.dx = 0.001
        self.dt = 0.001
        self.lattice_file = lattice_file
        self.hack_lattice()
        self.load_lattice()

    def hack_lattice(self):
        source = self.lattice_file
        src_dir = "/".join(self.lattice_file.split("/")[:-1])
        target_dir = src_dir+"/hack/"
        utils.utilities.clear_dir(target_dir)
        shutil.copy2(src_dir+"/disttest.dat", target_dir+"/disttest.dat")
        prefix = self.lattice_file.split(".")[:-1]+["hack"]
        target = ".".join(prefix)
        target = self.lattice_file.split("/")
        target.insert(-1, "hack")
        target = "/".join(target)
        print("Hacking from\n   ", source, "to\n   ", target)
        if source == target:
            raise RuntimeError("Source == target! oh noes!")
        fin = open(source, "r")
        fout = open(target, "w")
        for line in fin.readlines():
            if "ENABLEHDF5=" in line:
                line = "Option, ENABLEHDF5=False; // plotfields HACK!!! \n"
            if "REAL N_TURNS=" in line:
                line = "REAL N_TURNS=0.001; // plotfields HACK!!! \n"
            if "BOOL DO_MAGNET_FIELD_MAPS=" in line:
                line = "BOOL DO_MAGNET_FIELD_MAPS=False; // plotfields HACK!!! \n"
            fout.write(line)
        self.lattice_file = target

    def load_lattice(self):
        here = os.getcwd()
        a_dir, a_file = os.path.split(self.lattice_file)
        os.chdir(a_dir)
        PyOpal.parser.initialise_from_opal_file(a_file)
        os.chdir(here)

    def get_derivative(self, var1, var2, x, y, z, t):
        pos_vec = [x, y, z, t]
        var2 = ["x", "y", "z", "t"].index(var2)
        pos_vec[var2] += self.dx
        field_plus = self.get_field(*pos_vec)[var1]
        pos_vec[var2] -= 2*self.dx
        field_minus = self.get_field(*pos_vec)[var1]
        derivative = (field_plus-field_minus)/2.0/self.dx
        return derivative

    def get_div_b(self, x, y, z, t):
        div_b = self.get_derivative("bx", "x", x, y, z, t) + \
                self.get_derivative("by", "y", x, y, z, t) + \
                self.get_derivative("bz", "z", x, y, z, t)
        return div_b

    def get_curl_b(self, x, y, z, t):
        curl_b = [
            self.get_derivative("by", "z", x, y, z, t) - \
            self.get_derivative("bz", "y", x, y, z, t),
            self.get_derivative("bx", "z", x, y, z, t) - \
            self.get_derivative("bz", "x", x, y, z, t),
            self.get_derivative("bx", "y", x, y, z, t) - \
            self.get_derivative("by", "x", x, y, z, t)
        ]
        return curl_b

    def get(self, var, x, y, z, t):
        if var in self.field_vars:
            return self.get_field(x, y, z, t)[var]
        elif var == "div_b":
            return self.get_div_b(x, y, z, t)
        elif var == "curl_b":
            curl_b = self.get_curl_b(x, y, z, t)
            curl_b = (curl_b[0]**2+curl_b[1]**2+curl_b[2]**2)**0.5
            return curl_b


    def get_field(self, x, y, z, t):
        phi = math.atan2(y, x)
        oob, bx, by, bz, ex, ey, ez = \
                                PyOpal.field.get_field_value(x, y, z, t)
        btot = (bx**2+by**2+bz**2)**0.5
        br = bx*math.cos(phi)+by*math.sin(phi)
        bphi =  -bx*math.sin(phi)+by*math.cos(phi)

        etot = (ex**2+ey**2+ez**2)**0.5
        er = ex*math.cos(phi)+ey*math.sin(phi)
        ephi =  -ex*math.sin(phi)+ey*math.cos(phi)
        field = {"bx":bx, "by":by, "bz":bz, "btot":btot, "br":br, "bphi":bphi,
                 "ex":ex, "ey":ey, "ez":ez, "etot":etot, "er":er, "ephi":ephi,
                 "oob":oob}
        return field

    field_vars = ["bx", "by", "bz", "btot", "br", "bphi", "ex", "ey", "ez", "etot", "er", "ephi", "oob"]

class LoadOrbit(object):
    def __init__(self, file_name, allowed_id, test_function = None):
        self.file_name = file_name
        self.orbit = {}
        self.allowed_id = allowed_id
        self.test_function = test_function
        self.parse_file()
        self.r_phi_track_file()

    def parse_file(self):
        heading = self.heading
        fin = open(self.file_name)
        for item in heading:
            self.orbit[item] = []
        line = fin.readline().rstrip("\n")
        line = fin.readline().rstrip("\n")
        while line != "":
            line = fin.readline().rstrip("\n")
            words = line.split()
            if len(words) != len(heading):
                print("Line\n  "+line+"\nmismatched to heading\n  "+str(heading)+"\nin parse_file "+self.file_name)
            else:
                words = [self.types[i](x)*self.units[i] for i, x in enumerate(words)]
                if words[0] not in self.allowed_id:
                    continue
                is_okay = self.test_function == None or not self.test_function(words)
                if not is_okay:
                    print("Breaking due to failed test function at", words, self.test_function(words))
                    break
                for i, item in enumerate(heading):
                    self.orbit[item].append(words[i])
        print("Got", len(self.orbit["x"]), "lines from file "+self.file_name)

    @classmethod
    def fix_domain(self, phi):
        if phi < self.azimuthal_domain[0]:
            phi += (self.azimuthal_domain[1]-self.azimuthal_domain[0])
        elif phi >  self.azimuthal_domain[1]:
            phi -= (self.azimuthal_domain[1]-self.azimuthal_domain[0])
        return phi

    def r_phi_track_file(self):
        data = self.orbit
        data["r"] = list(range(len(data["x"])))
        data["phi"] = list(range(len(data["x"])))
        data["pr"] = list(range(len(data["x"])))
        data["pphi"] = list(range(len(data["x"])))
        for i in range(len(data["r"])):
            data["r"][i] = (data["x"][i]**2+data["y"][i]**2.)**0.5
            phi = math.atan2(data["y"][i], data["x"][i])
            data["phi"][i] = self.fix_domain(math.degrees(phi))
            px = data["px"][i]
            py = data["py"][i]
            data["pr"][i]   = px*math.cos(phi)+py*math.sin(phi)
            data["pphi"][i] = -px*math.sin(phi)+py*math.cos(phi)

    p_mass = xboa.common.pdg_pid_to_mass[2212]
    heading = ["id", "x", "px", "y", "py", "z", "pz"]
    units = [1, 1.0, p_mass, 1.0, p_mass, 1.0, p_mass]
    types = [str]+[float]*6
    azimuthal_domain = [-180.0, 180.0]

class Labels(object):
    labels = {
        "x":"x [m]", "y":"y [m]", "z":"z [m]",
        "r":"r [m]", "phi":"$\\phi$ [$^o$]",
        "px":"p$_{x}$ [MeV/c]", "py":"p$_{y}$ [MeV/c]", "pz":"p$_{z}$ [MeV/c]",
        "pr":"p$_{r}$ [MeV/c]", "pphi":"p$_{\\phi}$ [MeV/c]",
        "bx":"B$_{x}$ [T]", "by":"B$_{y}$ [T]", "bz":"B$_{z}$ [T]",
        "btot":"B$_{tot}$ [T]", "br":"B$_{r}$ [T]", "bphi":"B$_{\\phi}$ [T]",
        "ex":"E$_{x}$ [MV/m]", "ey":"E$_{y}$ [MV/m]", "ez":"E$_{z}$ [MV/m]",
        "etot":"E$_{tot}$ [MV/m]", "er":"E$_{r}$ [MV/m]", "ephi":"E$_{\\phi}$ [MV/m]",
        "div_b":"$\\nabla \\cdot \\mathbf{B}$ [T/m]",
        "curl_b":"$|\\nabla \\times \\mathbf{B}|$ [T/m]",
        "t":"Time [ns]"
    }

class RFAnalyser(object):
    def __init__(self, fields):
        self.fields = fields
        self.n_points = 100

    def get_data(self, centre, var, time_window):
        self.x_list = []

    def get_seed(self):
        pass

    def sine_fit(self, centre, var, time_window):
        voltage = max(self.y_list)
        frequency = 1e-3
        crossings = []
        for i, y in enumerate(self.y_list[1:]):
            if y > 0. and self.y_list[i] < 0: 
                crossings.append(i)
        if len(crossings) > 0:
            t0 = self.x_list[crossings[0]]
        print(crossings)#[1], crossings[0]
        #print "FIT CROSSINGS", self.x_list[crossings[1]], self.x_list[crossings[0]]
        if len(crossings) > 1:
            frequency = 1./(self.x_list[crossings[1]]-self.x_list[crossings[0]])
        frequency *= 2.*math.pi
        print("Seeding sine fit with", t0, frequency, voltage)
        fitter = ROOT.TF1("sin "+str(len(self.root_objects)), "[0]*sin([1]*(x-[2]))")
        fitter.SetParameter(0, voltage)
        fitter.SetParameter(1, frequency)
        fitter.SetParameter(2, t0)
        fitter.SetRange(min(self.x_list), max(self.x_list))
        #fitter.Draw("SAME")
        self.graph.Fit(fitter)
        self.canvas.Update()
        self.root_objects.append(fitter)
        rf_parameters = {
            "voltage":fitter.GetParameter(0),
            "frequency":fitter.GetParameter(1)/2./math.pi,
            "t0":fitter.GetParameter(2)
        }
        return rf_parameters

class PlotFields(object):
    def __init__(self, fields, orbit_plotter = []):
        self.n_1d_points = 1000
        self.n_2d_points = 1000
        self.fields = fields
        self.orbit_list = orbit_plotter
        self.cmap = "PiYG"
        self.r_min = 3.50
        self.r_max = 4.50
        self.n_phi = 10
        self.do_pipe = False

    def field_fig(self, job_name, centre, range_x, range_y, do_elevation=False):
        figure = matplotlib.pyplot.figure(figsize=(20, 10))

        axes = figure.add_subplot(2, 3, 1) 
        self.plot_2d(figure, axes, centre, range_x, range_y, [-1.0, 1.0], "x", "y", "bx")
        for orbit in self.orbit_list:
            orbit.plot_2d(axes, "x", "y", name=orbit.name)

        axes = figure.add_subplot(2, 3, 2) 
        self.plot_2d(figure, axes, centre, range_x, range_y, [-1.0, 1.0], "x", "y", "by")
        for orbit in self.orbit_list:
            orbit.plot_2d(axes, "x", "y", name=orbit.name)

        axes = figure.add_subplot(2, 3, 3) 
        self.plot_2d(figure, axes, centre, range_x, range_y, [-1.0, 1.0], "x", "y", "bz")
        for orbit in self.orbit_list:
            orbit.plot_2d(axes, "x", "y", name=orbit.name)

        axes = figure.add_subplot(2, 3, 4) 
        if not do_elevation:
            self.plot_1d(axes, centre, range_y, [-1.0, 1.0], "y", ["bx", "by", "bz", "btot", "div_b", "curl_b"], 0, "")
        else:
            for orbit in self.orbit_list:
                orbit.plot_phi_z(axes)

        axes = figure.add_subplot(2, 3, 5) 
        self.plot_2d(figure, axes, centre, range_x, range_y, [-0.1, 0.1], "x", "y", "div_b")
        for orbit in self.orbit_list:
            orbit.plot_2d(axes, "x", "y", name=orbit.name)

        axes = figure.add_subplot(2, 3, 6) 
        self.plot_2d(figure, axes, centre, range_x, range_y, [-0.1, 0.1], "x", "y", "curl_b")
        for orbit in self.orbit_list:
            orbit.plot_2d(axes, "x", "y", name=orbit.name)
        axes.legend()

        figure.suptitle(job_name)
        return figure


    def azimuthal_fig(self, job_name, range_phi):
        figure = matplotlib.pyplot.figure(figsize=(20, 10))
        axes = figure.add_subplot(2, 2, 1) 
        for orbit in self.orbit_list:
            orbit.plot_2d(axes, 'phi', 'z', range_phi, None)
        axes.legend()

        axes = figure.add_subplot(2, 2, 2) 
        for orbit in self.orbit_list:
            orbit.plot_2d(axes, 'phi', 'r', range_phi, None)
        axes.legend()

        axes = figure.add_subplot(2, 2, 3)
        self.plot_1d_orbit(axes, self.orbit_list[1], 'phi', ['br', 'bphi', 'bz', 'btot'], range_phi, None)

        figure.suptitle(job_name)
        return figure


    def rf_fig(self, job_name, centre, time_list, range_t):
        figure = matplotlib.pyplot.figure(figsize=(20, 10))
        figure.suptitle(job_name)
        axes = figure.add_subplot(1, 1, 1) 
        for t_offset in time_list:
            self.plot_1d(axes, centre, range_t, [-1.0, 1.0], "t", ["ephi", "er", "btot"], t_offset, "")

    def bump_fig(self, job_name, radius, phi_list, range_t, range_phi):
        figure = matplotlib.pyplot.figure(figsize=(20, 10))
        figure.suptitle(job_name)
        axes = figure.add_subplot(1, 2, 1) 
        for phi in phi_list:
            name = "phi: "+format(int(math.degrees(phi)))+"$^\\circ$"
            centre = [radius*math.cos(phi), radius*math.sin(phi), 0, 0]
            self.plot_1d(axes, centre, range_t, [-1.0, 1.0], "t", ["br", "bz"], 0, name)
        axes = figure.add_subplot(1, 2, 2) 
        for orbit in self.orbit_list:
            orbit.plot_phi_z(axes)
        axes.set_xlim(range_phi)
        axes.set_ylim([-0.23, -0.07])
        #self.traj_fig(axes)
        for orbit in self.orbit_list:
            orbit.plot_phi_z(axes)
        return figure

    def traj_fig(self, axes):
        min_z, max_z = -0.1, 0.1
        var_x = "phi"
        field_variables = ["br", "bz", "bphi"]
        phi_list = self.orbit_list[0].orbit.orbit["phi"]
        r_list = self.orbit_list[0].orbit.orbit["r"]
        min_elev, max_elev = axes.get_ylim()
        min_phi, max_phi = axes.get_xlim()
        step_phi = (max_phi-min_phi)/self.n_2d_points
        step_elev = (max_elev-min_elev)/self.n_2d_points
        track_index = 0
        x_list, y_list, z_list = [], [], []
        for i in range(self.n_2d_points-1):
            phi = step_phi*(i+0.5)+min_phi
            while phi < phi_list[track_index] and track_index < len(phi_list)-1:
                track_index += 1
            radius = r_list[track_index]
            x = math.cos(phi)*radius
            y = math.sin(phi)*radius
            for j in range(self.n_2d_points-1):
                elev = step_elev*(j+0.5)+min_elev
                z = self.fields.get("br", x, y, elev, 0.0)
                if z < min_z:
                    z = min_z
                if z > max_z:
                    z = max_z
                x_list.append(phi)
                y_list.append(elev)
                z_list.append(z)
        cmax = max(abs(min_z), abs(max_z))
        vtot = max([abs(min(z_list)), max(z_list)])
        vtot = cmax #min(vtot, cmax)
        hist = axes.hist2d(x_list, y_list, self.n_2d_points, [[min_phi, max_phi], [min_elev, max_elev],], False, z_list, 
                    cmin=min_z, cmax=max_z, cmap=self.cmap, vmin=-vtot, vmax=vtot)
        axes.get_figure().colorbar(hist[3], ax=axes)


    def plot_1d(self, axes, centre, range_x, range_y, var_x, var_y_list, offset_x, name):
        i_x = self.pos_vars.index(var_x)
        min_x, max_x = range_x[0], range_x[1]
        min_y, max_y = range_y[0], range_y[1]
        pos = copy.deepcopy(centre)
        step_x = (max_x-min_x)/self.n_1d_points
        x_list = [step_x*(i+0.5)+min_x for i in range(self.n_1d_points-1)]
        for var_y in var_y_list:
            y_list = []
            for i, x in enumerate(x_list):
                pos[i_x] = x+offset_x
                y_list.append(self.fields.get(var_y, pos[0], pos[1], pos[2], pos[3]))
            label = Labels.labels[var_y]
            if name:
                label += " "+name
            if offset_x:
                label += " offset "+format(offset_x, "6.4g")
            axes.plot(x_list, y_list, label=label)
            axes.set_xlabel(Labels.labels[var_x])
        axes.legend()

    def plot_2d(self, figure, axes, centre, range_x, range_y, range_z, var_x, var_y, var_z):
        i_x = self.pos_vars.index(var_x)
        i_y = self.pos_vars.index(var_y)
        min_x, max_x = range_x[0], range_x[1]
        min_y, max_y = range_y[0], range_y[1]
        min_z, max_z = range_z[0], range_z[1]
        pos = copy.deepcopy(centre)
        step_x = (max_x-min_x)/self.n_2d_points
        step_y = (max_y-min_y)/self.n_2d_points
        x_list = []
        y_list = []
        z_list = []
        for i in range(self.n_2d_points-1):
            x = step_x*(i+0.5)+min_x
            pos[i_x] = x
            for j in range(self.n_2d_points-1):
                y = step_y*(j+0.5)+min_y
                pos[i_y] = y
                z = self.fields.get(var_z, pos[0], pos[1], pos[2], pos[3])
                if z > max_z:
                    z = max_z
                if z < min_z:
                    z = min_z
                x_list.append(x)
                y_list.append(y)
                z_list.append(z)
        cmax = max(abs(min_z), abs(max_z))
        vtot = max([abs(min(z_list)), max(z_list)])
        vtot = cmax #min(vtot, cmax)
        hist = axes.hist2d(x_list, y_list, self.n_2d_points, [[min_x, max_x], [min_y, max_y]], False, z_list, 
                    cmin=min_z, cmax=max_z, cmap=self.cmap, vmin=-vtot, vmax=vtot)
        axes.set_xlabel(Labels.labels[var_x])
        axes.set_ylabel(Labels.labels[var_y])
        axes.set_title(Labels.labels[var_z])
        if self.do_pipe:
            self.cartesian_lines(axes)
        figure.colorbar(hist[3], ax=axes)

    def plot_1d_orbit(self, axes, orbit, var_x, var_y_list, range_x, range_y):
        i_list = [i for i, x in enumerate(orbit.orbit.orbit[var_x]) \
                                          if x < range_x[1] and x > range_x[0]]
        x_list = [orbit.orbit.orbit[var_x][i] for i in i_list]
        xcoord_list = orbit.orbit.orbit['x']
        ycoord_list = orbit.orbit.orbit['y']
        zcoord_list = orbit.orbit.orbit['z']
        y_dict_of_lists = dict([(var_y, []) for var_y in var_y_list])
        for i in i_list:
            x, y, z = xcoord_list[i], ycoord_list[i], zcoord_list[i]
            for var_y in var_y_list:
                y_dict_of_lists[var_y].append(self.fields.get(var_y, x, y, z, 0))
        for var_y, y_list in y_dict_of_lists.items():
            axes.plot(x_list, y_list, label=Labels.labels[var_y]) 
        axes.legend()
        axes.set_xlabel(Labels.labels[var_x])
        axes.set_title(orbit.name)
        if range_x != None:
            axes.set_xlim(range_x)
        if range_y != None:
            axes.set_ylim(range_y)

    def cartesian_lines(self, axes):
        for i_phi in range(self.n_phi):
            phi = math.pi*2.0/self.n_phi*i_phi
            x = self.r_min*math.cos(phi), self.r_max*math.cos(phi)
            y = self.r_min*math.sin(phi), self.r_max*math.sin(phi)
            axes.plot(x, y, color="gray")
        phi_list = [math.radians(i_phi) for i_phi in range(361)]
        for r0 in self.r_min, self.r_max:
            x = [r0*math.cos(phi) for phi in phi_list]
            y = [r0*math.sin(phi) for phi in phi_list]
            axes.plot(x, y, color="gray")

    pos_vars = ["x", "y", "z", "t"]

class LoadH5(object):
    def __init__(self, file_name_glob):
        self.time_window = 100
        self.height_window=[-0.2, 0.2]
        self.id_cut = [0]
        self.file_name_glob = file_name_glob
        self.data = []
        f_glob = glob.glob(self.file_name_glob)
        print ("Globbing ", self.file_name_glob, "yields", f_glob)
        for file_name in f_glob:
            self.load_h5_probe(file_name)

    def load_h5_probe(self, file_name):
        data = []
        print("Loading h5 file", file_name, end=" ")
        h5_file = h5py.File(file_name, 'r')
        for key in h5_file.keys():
            if key[:5] != "Step#":
                print("skipping key", key, end=" ")
                continue
            n_steps = len(h5_file[key]["x"])
            h5_step = h5_file[key]
            for i in range(n_steps):
                item = {
                    "id":h5_step["id"][i],
                    "x":h5_step["x"][i],
                    "y":h5_step["y"][i],
                    "z":h5_step["z"][i],
                    "t":h5_step["time"][i],
                    "phi":math.degrees(math.atan2(h5_step["y"][i],
                                                  h5_step["x"][i])),
                    "r":(h5_step["x"][i]**2+h5_step["y"][i]**2)**0.5
                }
                self.phi_domain(item)
                data.append(item)
                #print("H5 Loaded", item)
        if len(data) == 0 or len(data) == 1:
            return []
        data = data[1:]
        print("... found", len(data), "points and after cuts", end=" ")
        data = self.cut(data)
        print("...", len(data), "points")
        self.data += data

    @classmethod
    def phi_domain(cls, item):
        dphi = cls.azimuthal_domain[1]-cls.azimuthal_domain[0]
        while item["phi"] < cls.azimuthal_domain[0]:
                    item["phi"] += dphi
        while item["phi"] > cls.azimuthal_domain[1]:
                    item["phi"] -= dphi

    def cut(self, data):
        if self.id_cut:
            data = [item for item in data if item["id"] not in self.id_cut]
        if self.time_window and len(data):
            start_time = min([item["t"] for item in data])
            data = [item for item in data if item["t"] < start_time+self.time_window]
        if self.height_window:
            data = [item for item in data if item["z"] < self.height_window[1] and item["z"] > self.height_window[0]]
        return data

    azimuthal_domain = [-180.0, 180.0]

class PlotOrbit(object):
    def __init__(self, orbit, name="", probe_data=None):
        self.orbit = orbit
        self.probe = probe_data
        self.name = name

    def do_plots(self):
        self.plot_phi_z()

    def plot_phi_z(self, axes):
        self.plot_2d(axes, "phi", "z", "")

    def plot_probe(self, axes, var_x, var_y, range_x, range_y, color):
        x_list = [item[var_x] for item in self.probe.data]
        y_list = [item[var_y] for item in self.probe.data]
        #print("Plotting probes", x_list, y_list)
        axes.scatter(x_list, y_list, label=self.name, s=20, marker="o", facecolors="none", edgecolors=color)

    def plot_2d(self, axes, var_x, var_y, range_x, range_y):
        i_list = range(len(self.orbit.orbit[var_x]))
        i_list = [i for i in i_list \
            if self.orbit.orbit[var_x][i] < range_x[1] and self.orbit.orbit[var_x][i] > range_x[0]] 
        x_list = [self.orbit.orbit[var_x][i] for i in i_list]
        y_list = [self.orbit.orbit[var_y][i] for i in i_list]
        axes.set_xlabel(Labels.labels[var_x])
        axes.set_ylabel(Labels.labels[var_y])
        if range_x != None:
            axes.set_xlim(range_x)
        if range_y != None:
            axes.set_ylim(range_y)
        line2d = axes.plot(x_list, y_list, label=self.name)
        if self.probe != None:
            self.plot_probe(axes, var_x, var_y, range_x, range_y, color=line2d[0].get_color())

def bump_plot():
    ramp = "None"
    base_dir = "output/double_triplet_baseline/single_turn_injection/find_bump_parameters_20mm-vertical_v3=-0.08T_1"
    run_dir = os.path.join(base_dir, "tmp/track_bump")
    orbit_file_glob = os.path.join(run_dir, "VerticalSectorFFA-trackOrbit-*dat")
    lattice_file = os.path.join(run_dir, "VerticalSectorFFA.tmp")
    allowed_events = ["ID0", "ID1"]
    test_function = lambda words: words[5] > 0.1
    field = GetFields(lattice_file)
    orbit_plotter_list = []
    for orbit_file in glob.glob(orbit_file_glob):
        orbit = LoadOrbit(orbit_file, allowed_events, test_function)
        orbit_plotter_list.append(PlotOrbit(orbit))
    plot_fields = PlotFields(field, orbit_plotter_list)
    plot_fields.n_2d_points = 50
    range_x = [-5.0, 5.0]
    range_y = [-5.0, 5.0]
    job_name = "ramp "+ramp
    for z in [0.0]:
        centre = [0.0, 0.0, z, 0.0]
        figure = plot_fields.field_fig(job_name, centre, range_x, range_y, True)
        job_name = job_name.replace(" ", "_")
        figure.savefig(os.path.join(base_dir, job_name+"_z="+str(z)+".png"))

def plot_azimuthal_old():
    base_dir = "output/double_triplet_baseline/single_turn_injection/find_bump_parameters_x_20.0_y_0.0_mm_3"
    run_dir = os.path.join(base_dir, "tmp/track_bump")
    orbit_file_glob = os.path.join(run_dir, "VerticalSectorFFA-trackOrbit-*.dat")
    lattice_file = os.path.join(run_dir, "VerticalSectorFFA.tmp")
    allowed_events = ["ID0", "ID1"]
    test_function = lambda words: words[5] > 0.1
    orbit_list = []
    for orbit_file in glob.glob(orbit_file_glob):
        name = orbit_file.split("VerticalSectorFFA-trackOrbit")[1][:-4]
        orbit_list.append(PlotOrbit(LoadOrbit(orbit_file, allowed_events, test_function), name))
    plot_fields = PlotFields(GetFields(lattice_file), orbit_list)
    angle_range = [0.0, 360.0] #[104, 112] #  
    job_name = "closed_orbit_bump"
    figure = plot_fields.azimuthal_fig("Bump", angle_range)
    figure.savefig(os.path.join(base_dir, job_name+".png"))
    print("Files", orbit_file_glob, glob.glob(orbit_file_glob))

def plot_azimuthal():
    base_dir = "output/double_triplet_baseline/single_turn_injection/track_bump_parameters_x_0.0_y_30.0_mm_4"
    co_dir = "output/double_triplet_baseline/single_turn_injection/track_bump_parameters_x_0.0_y_20.0_mm_3.1/track_beam/forwards"
    run_dir = os.path.join(base_dir, "track_beam/")
    orbit_folder_list = [(os.path.join(run_dir, "forwards"), "bumped H$^{+}$"),
                         (os.path.join(run_dir, "backwards"), "injected H$^{-}$"),
                         (co_dir, "unbumped")]
    probe_files = "*PROBE*.h5" #"FOILPROBE_1.h5" #
    track_files = "VerticalSectorFFA-trackOrbit.dat"
    lattice_file = os.path.join(base_dir, "tmp/track_beam/VerticalSectorFFA.tmp")
    angle_domain = [-0.0, 360.0] #[104, 112] #  
    plot_range = [108-18, 108+18]
    allowed_events = ["ID1"]

    LoadOrbit.azimuthal_domain = angle_domain
    LoadH5.azimuthal_domain = angle_domain
    test_function = lambda words: words[5] > 0.1
    orbit_list = []
    for orbit_folder, name in orbit_folder_list:
        h5 = LoadH5(os.path.join(orbit_folder, probe_files))
        orbit_file = os.path.join(orbit_folder, track_files)
        #name = orbit_file.split("track_beam/")[1].split("/Vert")[0]
        orbit_list.append(PlotOrbit(LoadOrbit(orbit_file, allowed_events, test_function), name, h5))
    plot_fields = PlotFields(GetFields(lattice_file), orbit_list)
    job_name = "closed_orbit_bump"
    figure = plot_fields.azimuthal_fig("Bump", plot_range)
    figure.savefig(os.path.join(base_dir, job_name+".png"))

def plot_long_term():
    "Routine to plot things over a long tracking cycle"
    base_dir = "output/double_triplet_baseline/single_turn_injection/bump_scan/track_bump_r_10.0_theta_0.0_scan/"
    run_dir = os.path.join(base_dir, "track_beam/")
    probe_files = "RINGPROBE01.h5" #"FOILPROBE_1.h5" #
    track_files = "VerticalSectorFFA-trackOrbit.dat"
    lattice_file = os.path.join(base_dir, "tmp/track_beam/VerticalSectorFFA.tmp")
    angle_domain = [-0.0, 360.0] #[104, 112] #  
    plot_range = [0, 360]
    allowed_events = ["ID1"]
    orbit_folder_list = [(base_dir+"/track_beam/forwards_0", "0"), (base_dir+"/track_beam/forwards_0", "1"), ]

    LoadOrbit.azimuthal_domain = angle_domain
    LoadH5.azimuthal_domain = angle_domain
    test_function = lambda words: words[5] > 0.1
    orbit_list = []
    for orbit_folder, name in orbit_folder_list:
        h5 = LoadH5(os.path.join(orbit_folder, probe_files))
        orbit_file = os.path.join(orbit_folder, track_files)
        orbit_list.append(PlotOrbit(LoadOrbit(orbit_file, allowed_events, test_function), name, h5))
    plot_fields = PlotFields(GetFields(lattice_file), orbit_list)
    job_name = "da_scan"
    figure = plot_fields.azimuthal_fig("Bump", plot_range)
    figure.savefig(os.path.join(base_dir, job_name+".png"))

def bump_plot_t():
    base_dir = "output/double_triplet_baseline/single_turn_injection/track_beam_13-1"
    run_dir = os.path.join(base_dir, "tmp/find_closed_orbits")
    orbit_file = os.path.join(run_dir, "VerticalSectorFFA-trackOrbit.dat")
    lattice_file = os.path.join(run_dir, "VerticalSectorFFA.tmp")
    plot_fields = PlotFields(GetFields(lattice_file), None)
    job_name = "ramp"
    r0 = 3.74
    #centre = [r0*math.cos(6/10*2*math.pi), r0*math.sin(6/10*2*math.pi), 0.0, 0.0]
    figure = plot_fields.bump_fig("Bump", 3.74, [2*math.pi*i/10.0 for i in range(1, 6)], [-500, 5000])
    figure.savefig(os.path.join(base_dir, job_name+"_t.png"))

def rf_plot():
    base_dir = "output/double_triplet_baseline/single_turn_injection/field_test"
    run_dir = os.path.join(base_dir, "tmp/find_closed_orbits")
    orbit_file = os.path.join(run_dir, "VerticalSectorFFA-trackOrbit_1.dat")
    lattice_file = os.path.join(run_dir, "VerticalSectorFFA.tmp")
    plot_fields = PlotFields(GetFields(lattice_file), None)
    job_name = ""
    r0 = 3.74
    #centre = [r0*math.cos(6/10*2*math.pi), r0*math.sin(6/10*2*math.pi), 0.0, 0.0]
    figure = plot_fields.bump_fig("Bump", 3.74, [2*math.pi*i/10.0 for i in range(1, 6)], [0.0, 15000])

def main():
    plot_long_term()
    #bump_plot()

if __name__ == "__main__":
    main()
    matplotlib.pyplot.show(block = False)
    input("Press <CR> to finish")

