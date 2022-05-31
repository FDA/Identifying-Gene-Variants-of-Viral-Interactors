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

#Polymer is a graph model of proteins, where the nodes are amino acids, the edges are native contacts between amino acids, and the backbone contains edges between adjacent amino acids.

import networkx as nx
from itertools import combinations_with_replacement

class Polymer(nx.Graph):
	def __init__(self, list_of_edges, nresidues=None, backbone=None, distances=None):
		super().__init__(list_of_edges)
		if backbone == None:
			if nresidues == None:
				nodes = [u for u in sorted(self.nodes())]
			else:
				nodes = [u for u in range(nresidues)]
			self.backbone = set([])
			for i in range(len(nodes) - 1):
				self.backbone |= set([(nodes[i], nodes[i + 1])]) #set of tuples representing edges between adjacent residues
		else:
			self.backbone = backbone
		self._distances = distances

	@staticmethod
	def sorted_edge(edge):
		return tuple(sorted(edge))

	def sorted_edges(self):
		return [Polymer.sorted_edge(edge) for edge in self.edges()]

	def residues(self):
		return set(b[0] for b in self.backbone) | set(b[1] for b in self.backbone)

	def number_of_residues(self):
		return len(self.backbone) + 1

	def number_connected_components(self):
		return nx.number_connected_components(self)

	def edge_ordering(self):
		edges = [Polymer.sorted_edge(edge) for edge in self.edges()]
		return {edge : i for edge, i in zip(sorted(edges), range(len(edges)))}

	def distances(self):
		if self._distances == None:
			G = nx.Graph([edge for edge in self.edges()])
			for b in self.backbone:
				G.add_edge(*b)
			self._distances = {}
			shortest_paths = nx.shortest_path_length(G) #use shortest_path_length instead of shortest_path because we only need the length, and the paths themselves use a lot of memory
			shortest_paths = dict(shortest_paths)
			for u,v in combinations_with_replacement(G.nodes(), 2):
				self._distances[(u, v)] = self._distances[(v, u)] = \
					shortest_paths[u][v]
		return self._distances

	def segments(self):
		if 0 in self and len(self[0]) > 0:
			segments = [[0]]
		else:
			segments = []
		for u in range(1, self.number_of_residues()):
			if u in self and len(self[u]) > 0: #if u is a vertex and has neighbors
				if u - 1 not in self or len(self[u - 1]) == 0: #u-1 not a vertex or has no neighbors
					segments.append([]) #start new segment
				segments[-1].append(u)
		return segments

	@staticmethod
	def read(stream):
		G = nx.Graph()
		for line in stream:
			t = line.split()
			if len(t) == 0:
				continue
			if len(t) == 1:
				G.add_node(int(t[0]))
			else:
				try:
					G.add_edge(int(t[0]), int(t[1]))
				except (KeyError, TypeError):
					raise TypeError("ERROR: cannot interpret line in input structure file:%d" % line)
		backbone = set((i, i + 1) for i in range(max(G.nodes())))
		return Polymer([edge for edge in G.edges()], backbone=backbone)

	def write(self, stream):
		for node in self.residues():
			stream.write("%d\n" % node)
		for edge in self.edges():
			stream.write("%d %d\n" % edge)
