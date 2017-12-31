#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""Continous time random walk (ctrw), with fixed jump length, and time
step drawn from a power law distribution
"""

import numpy as np
import argparse

def parse_arg():
    parser = argparse.ArgumentParser(description="Continous time random walk")
    parser.add_argument('save_path', help="Path to save result to", type=str)
    parser.add_argument('-n', '--samplingpoints', help='Sampling points', required=False, type=int, default=100)
    parser.add_argument('-m', '--trajectories', help='Trajectories to use', required=True, type=int)
    parser.parse_args()
    return parser.parse_args()


def waiting_time():
    "Compute CTRW waiting time, power-law"

    r = 0
    while r == 0:
        r = np.random.rand()
    a = 1
    alpha = 0.5
    tau = a * (np.power(r, -1.0/alpha) - 1.0)
    return tau


def main(args):
    "Main program"

    M = int(args.trajectories)  # trajectories
    FUDGE = 10                  # space between time points to save
    N = int(args.samplingpoints)*FUDGE

    epsilon = 1                 # time step in simulation
    time = np.zeros((N,))       # all sampling times
    for i in range(N):
        time[i] = i * epsilon

    dx2 = np.zeros((M, N))      # save MSD, one dimension

    mu = 0                      # mean of gaussian jump distribution
    a = 1.0                     # lattice constant
    sigma = a

    for m in range(M):
        i = 0                   # index
        t = 0                   # current time
        pos = 0                 # current particle position

        while t < time[-1]:
            tau = waiting_time()

            while i < N and t <= time[i] and time[i] < t + tau:
                dx2[m, i] = pos**2
                i = i + 1

            pos = pos + a*np.random.normal(mu, sigma)
            #pos = pos + np.random.randint(0, 2)*2 -1
            t = t + tau

    mean = np.zeros((N,))
    for i in range(N):
        for m in range(M):
            mean[i] = mean[i] + dx2[m, i]
        mean[i] = mean[i] / M


    # for i in range(len(mean)):
    #     print "%s\t%s" % (time[i],mean[i])

    np.savez_compressed(args.save_path + "/trajectories", trajectories=dx2[:,::FUDGE], time=time[::FUDGE])
    return 0

if __name__ == "__main__":
    args = parse_arg()
    main(args)
