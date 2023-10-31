#!/usr/bin/env python
"""Initialize the project's data space.

Iterates over all defined state points and initializes
the associated job workspace directories.
The result of running this file is the creation of a signac workspace:
    - signac.rc file containing the project name
    - signac_statepoints.json summary for the entire workspace
    - workspace/ directory that contains a sub-directory of every individual statepoint
    - signac_statepoints.json within each individual statepoint sub-directory.

"""

import signac
import flow
import logging
from collections import OrderedDict
from itertools import product


def get_parameters():
    ''''''
    parameters = OrderedDict()
    parameters["forcefield"] = ["pps_opls"]
    parameters["num_mols"] = [10,20,30,40,50,60,70,80,90,100]
    parameters["lengths"] = [8]
    parameters["density"] = [1.35]
    parameters["remove_hydrogens"] = [True]
    parameters["remove_charges"] = [True]
    parameters["kT"] = [1.0]
    parameters["n_steps"] = [5e5]
    parameters["shrink_kT"] = [1.0]
    parameters["shrink_n_steps"] = [5e5]
    parameters["shrink_period"] = [10000]
    parameters["r_cut"] = [2.5]
    parameters["dt"] = [0.0003]
    parameters["shrink_tau_kT"] = [1.0]
    parameters["tau_kT"] = [0.01]
    parameters["gsd_write_freq"] = [10000]
    parameters["log_write_freq"] = [10000]
    return list(parameters.keys()), list(product(*parameters.values()))


def main():
    project = signac.init_project() # Set the signac project name
    param_names, param_combinations = get_parameters()
    # Create the generate jobs
    for params in param_combinations:
        statepoint = dict(zip(param_names, params))
        job = project.open_job(statepoint)
        job.init()
        job.doc.setdefault("sim_done", False)
        job.doc.setdefault("sample_done", False)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
