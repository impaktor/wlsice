#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""Oscillator with brownian noise, taken from
http://arxiv.org/abs/1102.0524
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


# np.float(Nx2) -> np.float(Nx2)
def simulate(X, time_step):
    """Function that does the main simulation. X is a Nx2 matrix,
    with elements X[0]=x, and X[1]=v"""
    N = len(X)

    kappa = 225e-6                # Hooke's constant, kg/s^2

    # critical
    #kappa = 10                   # Hooke's constant

    m = 1e-12                     # mass, kg
    kBT = 275 * 1.38064852e-23    # k_B * T

    # TEST setting constants to 1-ish:
    kappa = 1
    m = 1
    X[0,0] = 1
    kBT = 1e-2

    gamma_crit = np.sqrt(4 * m * kappa)    # critical dampening
    gamma = gamma_crit         # Friction cooefficient

    D = kBT / gamma               # Einstein relation
    dt = time_step                # Time step
    tau = m / gamma            # characteristic time of exp. decrease
    w_0 = np.sqrt(kappa/m)     # omega_0: cyclic freq. of un-damped oscilator

    if gamma == gamma_crit:
        # _really_ set omega to zero (or else numeric noise)
        w = 0
    else:
        w = np.sqrt(w_0**2-1.0/(4*tau**2)) # cyclic freq. of damped oscillator

    # Print true parameter values:
    # print "l1=", 1./(2*tau)
    # print "l2=", w

    # Keep w real
    # assert(gamma**2 <= gamma_crit)

    # At critical freq. w=0, and thus gamma_crit = gamma
    if w == 0 or gamma_crit == gamma:                  # At critical dampening
        sig2_xx = 4*D*tau* (1 - np.exp(-dt/tau)*(1 + dt/tau + 1./2 * (dt/tau)**2))
        sig2_vv = D/tau * (1 - np.exp(-dt/tau)*(1 - dt/tau + 1./2 * (dt/tau)**2))
        sig2_xv = D* np.exp(-dt/tau)*(dt/tau)**2
    else:
        koeff = D / (4*w**2*w_0**2*tau**3)
        sig2_xx = koeff *        (4*w**2*tau**2 + np.exp(-dt/tau)*(np.cos(2*w*dt)-2*w*tau*np.sin(2*w*dt)-4*w_0**2*tau**2))
        sig2_vv = koeff*w_0**2 * (4*w**2*tau**2 + np.exp(-dt/tau)*(np.cos(2*w*dt)+2*w*tau*np.sin(2*w*dt)-4*w_0**2*tau**2))
        sig2_xv = koeff*w_0**2*4*tau * np.exp(-dt/tau) * np.sin(w*dt)**2

    sig4_xv = sig2_xv**2
    sig_xx = np.sqrt(sig2_xx)
    sig_vv = np.sqrt(sig2_vv)
    sig_xv = np.sqrt(sig2_xv)

    # M =  np.array([[ 0.,  -1.], [ w_0**2,  1./tau]])

    mu = 0                      # mean of gaussian distribution
    sigma = 1                   # standard dev. of gaussian distribution

    # This is a constant, compute once:
    const = np.sqrt(sig2_vv - sig4_xv/sig2_xx)

    xi = np.random.normal(mu, sigma, N) # noise term
    ny = np.random.normal(mu, sigma, N) # noise term

    #print "# sig_xx, const:", sig_xx, const, xi, ny

    dX = np.zeros((N,2))
    dX[:,0] = sig_xx * xi
    dX[:,1] = np.add(sig2_xv/sig_xx * xi, const * ny)

    cosI = np.cos(w*dt) * np.array([[ 1.,  0.], [0.,  1.]])

    if w != 0 and gamma_crit != gamma:
        sinJ = np.sin(w*dt)*np.array([[ 1./(2*w*tau),  1./w], [ -w_0**2/w,  -1./(2*w*tau)]])
    else:
        sinJ = np.zeros((2,2))
        sinJ[0,0] =  dt/(2*tau)
        sinJ[0,1] = dt
        sinJ[1,0] = -w_0**2*dt
        sinJ[1,1] = - dt / (2*tau)

    eMdt = np.exp(-dt/(2*tau))* np.add(cosI, sinJ)      # 2x2 matix

    for i in range(len(X)-1):
        X[i+1] = np.add(np.dot(eMdt,X[i]), dX[i])
    return X



def main(args):
    N = args.samplingpoints
    M = args.trajectories

    traj = np.zeros((M, N))     # The data I want out from the program

    t_stop = 20.              # stop time is 20 ms.
    dt = t_stop / N             # time step size

    assert dt > 0

    A = 1e-9                    # Length we pull the spring at t_0, in meters
    for m in range(M):
        # This holds X = (x, v), i.e. position and speed
        X = np.zeros((N, 2))
        X[0,0] = A
        traj[m] = simulate(X, dt)[:,0] # save x-position

    time = np.zeros((N,))
    for i in range(N):
        time[i] = i * dt

    np.savez_compressed(args.save_path + "/trajectories", trajectories=traj, time=time)

    #Print mean
    mean = np.mean(traj,axis=0)
    print "#time\t mean"
    for i in range(len(mean)):
        print "%s\t%s" % (time[i],mean[i])

    return 0


if __name__ == "__main__":
    args = parse_arg()
    main(args)
