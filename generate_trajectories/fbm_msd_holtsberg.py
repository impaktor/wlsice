#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
#   --- Python-program fbm_msd_holtsberg.py ----
#
#   Simulate fractional Brownian motion
#   using an algorithm implmented by Holtsberg
#   [very similar to the Davies and Harte algorithm, see
#    R.B. Davies and D.S. Harte, Biometrika 74 (1987), p. 95;
#    M.J. Chambers, Math. Comput. Modelling 22, (1995), p. 1].
#
#   Input parameter H is 0 < H < 1, and it is 0.5 for
#   a standard Wiener process.
#
#   Note that we get two independent realisations for every
#   simulation.
#
#   This is version was ported to Pyton by Karl Fogelmark in 2014,
#   from the matlab version, by Tobias Ambjornsson, 2011, computing
#   the mean square displacement for FBM. This code is based on the
#   original code by Anders Holtsberg, 1999
#

import sys
import numpy as np
import matplotlib.pyplot as plt
import argparse

isPlottingOn = False
dumptrajectories = False        # Dump each trajectory as a text file

# ----- Input parameters -----
H = 0.25                   # Hurst exponent, 0<=H<1
# ----------------------------

def parse_arg():
    "Use argparse to parse the arguments from the command line"

    parser = argparse.ArgumentParser(description="Generate fractional Brownian motion data")
    parser.add_argument('-m','--trajectories', metavar="M", type=int, help='number of \"ensembles\"', required=True)
    parser.add_argument('-n','--sampling_points', metavar="N", type=int, help='closest (multiple of 2) number of sampling points', required=False, default = 2**7)
    parser.add_argument('-s','--stop_time', metavar="s", type=int, help='Stop time of simulation', required=False, default = 1000)
    parser.add_argument('-o','--output', help='output path', required=True, default = "out_")
    parser.add_argument('-b','--brute', help='use brute force', dest='bruteforce',
                        action="store_true", required=False)
    parser.set_defaults(bruteforce=False)
    args = parser.parse_args()
    return args


def generate_trajectories(no_ens, R, path, t, sampling_points, brute):
    N = sampling_points         # how many points to actually save
    mean_disp = np.zeros(N)     # mean displacement
    msd = np.zeros(N)           # mean square displaement
    n = len(t) - 2              # because reasons

    trajectories = np.zeros((no_ens, N))

    X = np.zeros(n)             # for brute force: store one point from each simulation run
    Y = np.zeros(n)             # for brute force: store one point from each simulation run
    for ens in range(no_ens):
        mu, sigma = 0, 1.0

        for i in range(n):      # re-run n times if brute force (else break out)

            # --- For H = 0 the process is white noise ---
            if H == 0:
                x = np.random.normal(mu,sigma,n)
                assert(len(x) == n)
                break

            # ------ Simulate  ---
            # Note the double length 2*n, that is the trick!
            x = np.random.normal(mu,sigma, 2*n) + complex(0,1)*np.random.normal(mu,sigma, 2*n)
            x = np.fft.ifft(np.multiply(np.sqrt(R), x)) * np.sqrt(2*n)
            x = x[0:n]

            # ------- Sum it up ---
            y = np.cumsum(np.imag(x))
            x = np.cumsum(np.real(x))

            if brute:
                X[i] = x[i]
                Y[i] = y[i]
            else:
                X = x
                Y = y
                break

        # --Remove redundant time points --
        T = np.zeros(N)
        for i in range(N):
            X[i] = X[int(i * n/N)]
            Y[i] = Y[int(i * n/N)]
            T[i] = t[int(i * n/N)]
        x = np.resize(X, N)
        y = np.resize(Y, N)

        # -- Mean and Mean Square displacement -------
        mean_disp = mean_disp + x + y
        msd = msd + np.power(x, 2) + np.power(y, 2)

        trajectories[ens] = np.add(np.power(x, 2), np.power(y, 2))

        if dumptrajectories:
            trajectory = np.add(np.power(x, 2), np.power(y, 2))
            name = path + ("/out_%s" % ens)
            f = open(name, 'w')
            for j in range(len(trajectory)):
                f.write("%s\t%s\n" % (T[j], trajectory[j]))
            f.close()

    # save all trajectories as a single compressed binary file
    np.savez_compressed(path + "/trajectories", trajectories=trajectories, time=T)

    # -- Normalize --
    mean_disp = mean_disp / (2 * no_ens)
    msd = msd / (2 * no_ens)
    T = np.zeros(N)
    for i in range(N):
        T[i] = t[int(i * n/N)]
    return (msd, mean_disp, T)


def plot(t, n, msd, msd_disp):
    a = 2*H
    msd_theory = np.power(t[0:n],a)

    if isPlottingOn:
        plt.loglog(t[0:n],msd, 'b', label="fBm")
        #plt.loglog(t[0:n],mean_disp, 'r')
        plt.loglog(t[0:n],msd_theory, 'g--', label="theory")
        plt.xlabel('t')
        plt.ylabel('MSD')
        plt.savefig('fbm.png', bbox_inches='tight')


def main(args):
    # -- convert input to easier to handle for FFT --
    n = 2**round(np.log(args.stop_time)/np.log(2))
    t = np.linspace(0,n+1,n+2)     # [0,1,2,....n,n+1], in increments of 1
    assert(len(t) == n + 2)          # since added a 0 in the beginning and a 1 in the end
    assert(args.sampling_points <= args.stop_time)

    # -- Compute the autocovariance function ---
    a = 2*H
    r = np.subtract((np.power(t+1,a) + np.power(np.abs(t-1),a))/2, np.power(t,a))  # TODO: probably should be /2.0

    # ----- The spectrum R ---
    r = np.append(r, r[n-2:0:-1])
    R = np.real(np.fft.fft(r))

    # ---- Check nonnegative R ---
    if np.any(R < 0):
        sys.stderr.write('This can happen only due to roundoff errors\n')

    # -- Ensemble-average loop  --
    mean_disp = np.zeros(n)     # mean displacement
    msd = np.zeros(n)           # mean square displaement
    (msd, mean_disp, T) = generate_trajectories(args.trajectories, R, args.output, t, args.sampling_points, args.bruteforce)

    # -------- Save ----------
    name = args.output + "/msd"
    f = open(name, 'w')
    for i in range(len(msd)):
        f.write("%s\t%s\n" % (T[i], msd[i]))
    f.close()

    # -------- Plot ----------
    if isPlottingOn:
        plot(T, n, msd, msd_disp)


if __name__ == "__main__":
    args = parse_arg()
    main(args)
