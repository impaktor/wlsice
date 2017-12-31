#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""Brownian motion, with fixed time step (epsilon), and jump lenght
drawn from gaussian distribution"""

import numpy as np
import argparse

def parse_arg():
    parser = argparse.ArgumentParser(description="Brownian motion random walk")
    parser.add_argument('save_path', help="Path to save result to", type=str)
    parser.add_argument('-n', '--samplingpoints', help='Sampling points', required=False, type=int, default=100)
    parser.add_argument('-m', '--trajectories', help='Trajectories to use', required=True, type=int)
    parser.parse_args()
    return parser.parse_args()


def main(args):
    "Main program"

    epsilon = 1                 # timestep
    mu = 0                      # mean of gaussian jump distribution
    a = 1.0                     # lattice constant
    sigma = a
    M = int(args.trajectories)
    FUDGE = 10
    N = int(args.samplingpoints)*FUDGE

    dx2 = np.zeros((M, N))
    for m in range(M):
        steps = a*np.random.normal(mu, sigma, N)
        # uniform_steps = np.random.randint(0, 2, N)*2 -1 # -1 or 1
        # steps = uniform_steps
        steps[0] = 0            # first step is in origin, at t=0
        dx2[m] = np.square(np.cumsum(steps))

        # f = open(args.save_path + "/out_%s" % m, 'w')
        # for i in range(len(dx2[m])):
        #     f.write("%s\t%s\n" % (i*epsilon, dx2[m,i]))
        # f.close()

    time = np.zeros((N,))
    for i in range(N):
        time[i] = i * epsilon

    np.savez_compressed(args.save_path + "/trajectories", trajectories=dx2[:,::FUDGE], time=time[::FUDGE])

    return 0

if __name__ == "__main__":
    args = parse_arg()
    main(args)
