import os
import shutil
import sys
import importlib

import ROOT

import xboa.common

import analysis.find_closed_orbits_4d
import analysis.find_tune
import analysis.find_da
import analysis.find_bump_parameters
import analysis.track_bump

from utils import utilities

def get_config():
    if len(sys.argv) < 2:
        print("Usage: python /path/to/run_one.py /path/to/config.py")
        sys.exit(1)
    config_file_name = sys.argv[1]
    config_file = sys.argv[1].replace(".py", "")
    config_file = config_file.split("scripts/")[1]
    config_file = config_file.replace("/", ".")
    print("Using configuration module", config_file)
    config_mod = importlib.import_module(config_file)
    config = config_mod.Config()
    return config_file_name, config

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

def master_substitutions(config):
    xboa.common.substitute(config.tracking["master_file"], config.tracking["lattice_file"], config.master_substitutions)

def main():
    config_file_name, config = get_config()
    output_dir(config, config_file_name)
    utilities.setup_gstyle()
    ROOT.gErrorIgnoreLevel = config.run_control["root_verbose"]
    #master_substitutions(config)
    if config.run_control["find_closed_orbits_4d"]:
        analysis.find_closed_orbits_4d.main(config)
    if config.run_control["find_tune"]:
        analysis.find_tune.main(config)
    if config.run_control["find_da"]:
        analysis.find_da.main(config)
    if config.run_control["find_bump_parameters"]:
        analysis.find_bump_parameters.main(config)
    if config.run_control["track_bump"]:
        analysis.track_bump.main(config)

if __name__ == "__main__":
    main()
