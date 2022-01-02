import os
import shutil
import sys
import importlib
import subprocess
import datetime

import numpy
import ROOT

import xboa.common

import analysis.find_closed_orbits_4d
import analysis.find_tune
import analysis.find_da
import analysis.find_bump_parameters
import analysis.build_bump_surrogate_model
import analysis.track_bump
import analysis.track_beam

from utils import utilities

def output_dir(config, config_file_name):
    output_dir = config.run_control["output_dir"]
    if config.run_control["clean_output_dir"]:
        try:
            shutil.rmtree(output_dir)
        except OSError:
            pass
    try:
        os.makedirs(output_dir)
    except OSError:
        pass
    shutil.copy2(config_file_name, output_dir)
    git_string = "Time "+str(datetime.datetime.now())+"\n"
    try:
        git_string += subprocess.check_output(["git", "log", "-1"]).decode('unicode_escape')
        git_string += subprocess.check_output(["git", "status"]).decode('unicode_escape')
    except Exception:
        git_string += "Error calling git"
    fout = open(output_dir+"/status", "w")
    print(git_string, file=fout)

def master_substitutions(config):
    xboa.common.substitute(config.tracking["master_file"], config.tracking["lattice_file"], config.master_substitutions)

def main():
    config_file_name, config = utilities.get_config()
    output_dir(config, config_file_name)
    utilities.setup_gstyle()
    ROOT.gErrorIgnoreLevel = config.run_control["root_verbose"]
    if config.run_control["random_seed"] != None:
        numpy.random.seed(config.run_control["random_seed"])
    #master_substitutions(config)
    if config.run_control["find_closed_orbits_4d"]:
        analysis.find_closed_orbits_4d.main(config)
    if config.run_control["find_tune"]:
        analysis.find_tune.main(config)
    if config.run_control["find_da"]:
        analysis.find_da.main(config)
    if config.run_control["find_bump_parameters"]:
        analysis.find_bump_parameters.main(config)
    if config.run_control["build_bump_surrogate_model"]:
        analysis.build_bump_surrogate_model.main(config)
    if config.run_control["track_bump"]:
        analysis.track_bump.main(config)
    if config.run_control["track_beam"]:
        analysis.track_beam.main(config)
    print("Finished with output in", config.run_control["output_dir"])

if __name__ == "__main__":
    try:
        main()
    except:
        raise
    finally:
        print('\033[0m')
              