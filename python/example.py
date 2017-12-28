#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Example file for using WLS-ICE method to fit a non-linear /
powerlaw function to mean ensamble data.

"""

import sys
import os
import numpy as np
import glob                     # For reading in data
import wlsice

# np.array(1), np.array(1) -> np.array(1)
def f(t, params):
    """A powerlaw: a*t**b"""
    assert(len(params) == 2)
    return np.multiply(params[0], np.power(t, params[1]))


# np.array(1), np.array(1) -> np.array(2)
def df(t, params):
    """Gradient of a powerlaw, with regards to parameters a,b:
    gradient(a*t**b) = [t**b, a*log(t)*t**b]
    """
    assert(t[0] > 0)
    assert(len(params) == 2)
    a = params[0]
    b = params[1]

    dfd_lam = np.zeros((len(params), len(t)))
    dfd_lam[0] = np.power(t, b) # list of N numbers
    dfd_lam[1] = np.multiply(a*np.log(t), np.power(t, b)) # list of N numbers
    return dfd_lam

# np.array(1), np.array(1) -> np.array(3!)
def d2f(t, params):
    """Return Hessian of powerlaw"""
    a = params[0]
    b = params[1]

    d2f = np.zeros((len(params), len(params), len(t)))

    # d2f / da da
    # df2[0,0] = np.zeroes(np.shape(t)) #(already there)

    # d2f / da db
    d2f[0,1] = np.multiply(np.log(t), np.power(t,b))

    # d2f / db da
    d2f[1,0] = d2f[0,1]

    # d2f / db db
    d2f[1,1] = np.multiply(a*(np.log(t)**2), np.power(t,b))

    return d2f


# string -> nil
def error(string, stop=True):
    sys.stderr.write("- Error! " + string)
    if stop:
        sys.exit(2)

# string, string -> tuple of [np.array(1), np.array(2)]
def read_in_data(load_base_path):
    """Read in all M 2-column data files (trajectories) from path,
    return a list of N time points and an MxN matrix of corresponding
    trajectory values.

    (Code that is commented out allows for saving the data back as
    a single compressed binary file, useful if many load/read operations)
    """

    trajectories = []
    time = []

    if not load_base_path:
        error("Specify base path to folder with data to fit. Each file having two colums, <time, trajectory>.\n")

    fileNames = glob.glob(load_base_path + "/*")
    # fileNames = os.listdir("load_base_path")

    if not fileNames:
        error("Specify <base path> to data files. No files found!\n")

    for i, fileName in enumerate(fileNames):
        f = open(fileName,'r')
        trajectories.append([])
        for line in f.readlines():
            li = line.strip()               # remove leading whitespace
            if not li.startswith("#"):      # ignore comments
                tmp = [float(value) for value in line.split()]
                trajectories[i].append(tmp[1])
                if i == 0:                  # save the time first time around
                    time.append(tmp[0])
        f.close()

    return np.array(time), np.array(trajectories)


# np.array(1), np.array(1) -> nil
def coff_det(y, f):
    """Coefficient of determination, to check goodness-of-fit, from mean
    of data (MSD), and f(t, params).
    https://en.wikipedia.org/wiki/Coefficient_of_determination
    """
    N = len(y)

    # note, by y_mean here we mean the average over N, not M!
    y_mean = np.sum(y) / N

    SS_tot = np.sum(np.square(np.subtract(y_mean, y)))
    SS_reg = np.sum(np.square(np.subtract(y_mean, f)))
    SS_res = np.sum(np.square(np.subtract(y,      f)))

    print("# Coefficient of determination:\t%s" % (1- SS_res / SS_tot))



def main(args):
    "Perform a least squares including correlation in error on data in folder"

    if len(args) < 2:
        print("Usage: %s <data-path/>" % args[0])
        sys.exit(1)

    # Read in data. Assumes path given is a folder where _all_ files
    # are 2 column files with x ("time") being column 1 and y
    # ("trajectory") being column 2..
    time, trajectories = read_in_data(args[1])
    M, N = np.shape(trajectories)

    min_method = "nm"           # Use Nelder-Mead minimization method
    guess=np.array([.5,.5])     # Starting point in parameter space

    #### OPTIONAL: Find index of first time to include in fitting procedure
    starttime = 200              # first time point to include in fitting
    epsilon0 = time[1] - time[0] # time step size
    start_idx = int(starttime / epsilon0)

    trajectories = trajectories[:,start_idx:]
    time = time[start_idx:]
    M,N = np.shape(trajectories)
    assert(N == len(time))
    print("# Removed the first %s sampling times, new N = %s" % (start_idx, N))
    #### END


    #### OPTIONAL: Reduce M (if we want to work with a smaller set of trajectories)
    new_M = 0                   # If 0, don't use this
    if new_M > 0:
        print("# Found M = %s trajectories" % M)
        if M < new_M:
            error("Have: M = %s, requested: %s\n" % (M, new_M))
        start_traj = 0          # can be changed to use a _different_ set of the full M
        stop_traj = start_traj + new_M
        assert(stop_traj <= M)
        trajectories = trajectories[start_traj:stop_traj, :]
        M = new_M
    #### END


    # The analytical function to fit, its gradient, and hessian
    wlsice.init(f, df, d2f)

    # Perform the actual fit
    params, sigma, chi2_min = wlsice.fit(time, trajectories, guess, min_method)

    # RESULT:
    print("# trajectories M=%s,\tsampling times N=%s, t_0=%s" % (M, N, time[0]))
    print("# Optimal param:\t%s" % params)
    print("# Sigma:\t%s" % sigma)
    print("# Chi-square value: \t%s" % chi2_min)

    # Also get goodness-of-fitt parametes
    y_mean = wlsice.computeMean(trajectories)
    coff_det(y_mean, f(time, params))


if __name__ == "__main__":
    main(sys.argv)
