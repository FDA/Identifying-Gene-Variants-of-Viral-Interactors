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
from itertools import combinations, product

from polymer import Polymer

class SubstructureError(Exception):
    pass

#returns connected components of line/dual graph, which frequently correspond to secondary structures
def find_substructures(polymer, min_component_size_conditional=6, min_component_size=10, nbr_dist=1):
    G = nx.Graph()
    for edge1 in polymer.edges():
        G.add_node(edge1) #every native contact is a node
        for edge2 in polymer.edges():
            if edge1 >= edge2: continue #edge2 hasn't been added to G yet, so don't try to add an edge involving it
            edgea = sorted(edge1)
            edgeb = sorted(edge2)
            if math.fabs(edgea[0] - edgeb[0]) <= nbr_dist and \
               math.fabs(edgea[1] - edgeb[1]) <= nbr_dist: #vertices in edges are close
                G.add_edge(edge1, edge2) #edge between native contacts involving neighboring amino acids.
    substructures = [cc for cc in reversed(sorted(nx.connected_components(G), key=len)) \
                     if (len(cc) >= min_component_size_conditional \
                         and len(set(v for edge in cc for v in edge)) <= len(cc)) \
                     or (len(cc) >= min_component_size)] #connected components of G of sufficient size
    component_edges = sorted([list(set(s for s in cc)) for cc in substructures], \
                             key=lambda y: min(x for x in y)) #list of (list of vertices) for each connected component
    return component_edges

def calc_ghost_edges(polymer):
    substructures = find_substructures(polymer)
    ghost_edges = []
    def attempt_add_edge(G, u, v): #add edges between u and neighbors of v, if u neighbors a vertex next to the neighbor on the backbone
        for nbr in G[v]:
            for nbrnbr in [nbr - 1, nbr + 1]:
                if nbrnbr in G[u]:
                    G.add_edge(nbr, u)
                    ghost_edges.append(tuple(sorted((nbr, u))))
                    return
    for substructure in substructures:
        G = nx.Graph([edge for edge in substructure])
        if nx.number_connected_components(G) > 1:
            for u in sorted(G.nodes()):
                ccs = [cc for cc in nx.connected_components(G)]
                cc = {u : c for c in range(len(ccs)) for u in ccs[c]}
                if u + 1 in G and cc[u + 1] != cc[u]:
                    attempt_add_edge(G, u, u + 1)
    G = nx.Graph(polymer.edges())
    if nx.number_connected_components(G) > 1:
        for u in sorted(G.nodes()):
            ccs = [cc for cc in nx.connected_components(G)]
            cc = {u : c for c in range(len(ccs)) for u in ccs[c]}
            if u + 1 in G and cc[u + 1] != cc[u]:
                attempt_add_edge(G, u, u + 1)
    return ghost_edges

#computes all subsets of substructures where all substructures are connected by edges
def calc_independent_configs(polymer, substructures):
    component_indices = [i for i in range(len(substructures))]
    independent_configs = []
    component_interactions = []
    unassigned_edges = [tuple(sorted(edge)) for edge in polymer.edges() \
                        if not any(edge in substructure for substructure in substructures)] #edges in the polymer not included in any substructure
	#find all edges between different substructures
    for i in component_indices:
        for j in component_indices:
            if i >= j: continue
            for edge1, edge2 in product(substructures[i], substructures[j]):
                if any((u == v or tuple(sorted((u,v))) in unassigned_edges) \
                       for u,v in [(edge1[0], edge2[0]), (edge1[1], edge2[0]), \
                                   (edge1[0], edge2[1]), (edge1[1], edge2[1])]):
                    component_interactions.append((i,j))
                    component_interactions.append((j,i))
                    break
	#generate all subsets of substructures such that all substructures are connected by edges
    for n in range(1, len(substructures) + 1):
        for group in combinations(component_indices, n): 
            G = nx.Graph([(i,j) for i,j in combinations(group, 2) if (i,j) in component_interactions]) #graph consisting of all substructures in group connected by edges
            for i in group:
                G.add_node(i)
            if nx.number_connected_components(G) == 1: #all substructures in group are connected by component_interactions
                independent_configs.append(set(group))
    return independent_configs

def get_polymer_segments(polymer, substructures):
    return {c : Polymer(substructures[c], backbone=polymer.backbone).segments() \
            for c in range(len(substructures))}

def get_expected_segments(polymer, polymer_segments, substructures, config):
    expected = {c : polymer_segments[c] for c in config}
    for c,segs in expected.items():
        expected[c] = [[u for u in seg if u in polymer and \
                        any(tuple(sorted((u, v))) in substructures[c] for v in polymer[u])] \
                       for seg in segs]
    return {c : segs for c,segs in expected.items() if all(len(seg) > 0 for seg in segs)}

def get_config_vertices(polymer, substructures, config):
    ss_vertices = set(u for c in config for edge in substructures[c] for u in edge)
    if len(ss_vertices) == 0:
        return set()
    min_ss_u, max_ss_u = min(ss_vertices), max(ss_vertices)
    other_vertices = set(u for u in polymer.nodes() if u > min_ss_u and u < max_ss_u) \
                     - set(u for c in range(len(substructures)) for edge in substructures[c] \
                           for u in edge)
    all_vertices = ss_vertices | other_vertices
    return all_vertices

def get_polymer_from_config(polymer, polymer_segments, substructures, config, kuhn_length, DEBUG=False):
    ss_vertices = set(u for c in config for edge in substructures[c] for u in edge)
    if len(ss_vertices) == 0:
        return Polymer([], nresidues=polymer.number_of_residues())
    min_ss_u, max_ss_u = min(ss_vertices), max(ss_vertices)
    other_vertices = set(u for u in polymer.nodes() if u > min_ss_u and u < max_ss_u) \
                     - set(u for c in range(len(substructures)) for edge in substructures[c] \
                           for u in edge)
    all_vertices = ss_vertices | other_vertices
    G = nx.Graph([edge for edge in polymer.edges() \
                  if edge[0] in all_vertices and edge[1] in all_vertices])
    cc_vertices = [u for cc in nx.connected_components(G) if any(v in cc for v in ss_vertices) \
                   for u in cc]
    subpolymer = Polymer([edge for edge in polymer.edges() \
                          if edge[0] in cc_vertices and edge[1] in cc_vertices], \
                         nresidues=polymer.number_of_residues())
    return subpolymer

def calc_co(component):
    return sum(math.fabs(edge[1] - edge[0]) for edge in component) / len(component)

def calc_fraction_helix(component):
    return sum(1 if math.fabs(edge[1] - edge[0]) < 4.5 else 0 for edge in component) / len(component)
