# Copyright (C) 2017 William M. Jacobs

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

import math
import networkx as nx
import numpy as np
import scipy.optimize

#find all sets of vertices not contained in segments
def find_loops(polymer, segments):
    connected_components = [cc for cc in nx.connected_components(polymer)]
    vertices = set(u for segs in segments.values() for seg in segs for u in seg)
    if len(vertices) > 0:
        minu, maxu = min(vertices), max(vertices)
        loops = [[]]
        for u in range(minu + 1, maxu + 1):
            if u - 1 not in vertices and u in vertices:
                loops.append([])
            elif u not in vertices and len(loops) > 0:
                loops[-1].append(u)
        segment_loops = [loop for loop in loops[:-1] if any(loop[0] - 1 in cc and loop[-1] + 1 in cc \
                                                            for cc in connected_components)]
    else:
        segment_loops = []
    return segment_loops

#computes S11
def loop_entropy(u, v, distances, mu, d, kuhnlength):
    l = v - u - 1
    if l <= kuhnlength:
        return -mu * l
    else:
        if distances[(round(u),round(v))] == 1:
            r = 0
        else:
            r = distances[(round(u),round(v))] / kuhnlength
        return -mu * kuhnlength + \
            d / 2. * (math.log(l / kuhnlength) + r**2 / l)

#returns minimum free energy for polymer
def minimize_fe(polymer, epsilon, distances, segments, mu=-2, d=3, kuhnlength=2, \
                l1=-0.8, l2=-0.4, g=0, nstar=0, verbose=False):
    warning = False
    n = epsilon.shape[0]
    zeronodes = set(u for u in range(n) if not np.count_nonzero(epsilon[u,:]))
    nodemask = np.array([1 if u not in zeronodes else 0 for u in range(n)])
    loops = find_loops(polymer, segments)
    t_epsilon = np.copy(epsilon) #S16
    t_mu = np.zeros(n) #S17
    for u in range(n):
        if u not in zeronodes:
            t_mu[u] = mu - 2 * g * nstar
            for v in range(n):
                if v not in zeronodes:
                    t_epsilon[u,v] += -g
                    t_epsilon[v,u] += -g
                    if v == u - 1 or v == u + 1:
                        t_epsilon[u,v] += l1
                        t_epsilon[v,u] += l1
                        t_mu[u] += l1
                    if (v == u - 2 and u - 1 not in zeronodes) or \
                       (v == u + 2 and u + 1 not in zeronodes):
                        t_epsilon[u,v] += l2
                        t_epsilon[v,u] += l2
                        t_mu[u] += l2
    def boltz(rho): #computes S15
        return np.exp(t_mu - np.dot(t_epsilon, rho)) * nodemask
    def dF_diff(rho):
        b = boltz(rho)
        return (rho - b / (1 + b)) * nodemask
    def free_energy(rho, nloop=4): #computes S14
        b = boltz(rho)
        F_contact = -0.5 * np.dot(rho, np.dot(t_epsilon, rho)) - np.dot(np.log(1 + b), nodemask)
        F_loop = 0
        for loop in loops:
            F_loop_i = 0
            norm_i = 0
            w_u = 1.
            for u in range(loop[0] - 1, max(-1, loop[0] - nloop), -1):
                w_u *= 1 - rho[u + 1]
                w_v = 1.
                for v in range(loop[-1] + 1, min(n, loop[-1] + nloop)):
                    w_v *= 1 - rho[v - 1]
                    dS = loop_entropy(u, v, distances, mu, d, kuhnlength)
                    w = w_u * w_v * rho[u] * rho[v] * math.exp(-dS)
                    F_loop_i += w * dS
                    norm_i += w
            F_loop += F_loop_i / norm_i
        return F_contact + F_loop - g * nstar**2
    if g == 0:
        rho_guess = np.ones(n)
    else:
        rho_guess = nstar / (n - len(zeronodes)) * np.ones(n)
    for u in zeronodes:
        rho_guess[u] = 0
    result = scipy.optimize.root(dF_diff, rho_guess, method='hybr',
                                 options={'xtol' : 1.e-12, 'maxfev' : 10000})
    if not np.allclose(dF_diff(result.x), np.zeros(n), atol=1.e-4):
        warning = True
        if verbose: print("WARNING: ftol = %g" % max(dF_diff(result.x)))
    if any(x < 0 or x > 1 for x in result.x):
        warning = True
        if verbose: print("WARNING: ! 0 < rho < 1")
    return result.x, free_energy(result.x), warning

# computes probability of contact formation between residues u and v
def avg_contacts(epsilon, rho):
    n = epsilon.shape[0]
    C = np.zeros((n, n))
    for u in range(n):
        for v in range(n):
            if epsilon[u,v] != 0:
                C[u,v] = rho[u] * rho[v]
    return C


def segment_end_probability(epsilon, rho, segments):
    avg_seg_end_prob = {}
    for c in segments:
        for i in range(len(segments[c])):
            for end,other,sign in [(0, -1, 1), (-1, 0, -1)]:
                w_v = 1.
                end_v = 0.
                norm_v = 0.
                for v in range(segments[c][i][end], segments[c][i][other] + sign, sign):
                    if v - sign >= 0 and v - sign < len(rho):
                        w_v *= 1 - rho[v - sign]
                    w = w_v * rho[v]
                    end_v += v * w
                    norm_v += w
                avg_seg_end_prob[(c, i, end)] = end_v / norm_v
    return avg_seg_end_prob

def subpolymer_epsilon(epsilon, vertices):
    epsilon_config = np.zeros(epsilon.shape)
    for u in vertices:
        for v in vertices:
            epsilon_config[u,v] = epsilon[u,v]
    return epsilon_config
