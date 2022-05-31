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

residue_abbreviations = \
{
    'G' : 'gly', # nonpolar
    'A' : 'ala',
    'V' : 'val',
    'L' : 'leu',
    'I' : 'ile',
    'M' : 'met',
    'F' : 'phe',
    'W' : 'trp',
    'P' : 'pro',
    'S' : 'ser', # polar
    'T' : 'thr',
    'C' : 'cys',
    'Y' : 'tyr',
    'N' : 'asn',
    'Q' : 'gln',
    'D' : 'asp', # negative
    'E' : 'glu',
    'K' : 'lys', # positive
    'R' : 'arg',
    'H' : 'his'
}
residue_abbreviations_inv = {v : k for k,v in residue_abbreviations.items()}

codon_table = {'I' : ['ATT', 'ATC', 'ATA'],
               'L' : ['CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG'],
               'V' : ['GTT', 'GTC', 'GTA', 'GTG'],
               'F' : ['TTT', 'TTC'],
               'M' : ['ATG'],
               'C' : ['TGT', 'TGC'],
               'A' : ['GCT', 'GCC', 'GCA', 'GCG'],
               'G' : ['GGT', 'GGC', 'GGA', 'GGG'],
               'P' : ['CCT', 'CCC', 'CCA', 'CCG'],
               'T' : ['ACT', 'ACC', 'ACA', 'ACG'],
               'S' : ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
               'Y' : ['TAT', 'TAC'],
               'W' : ['TGG'],
               'Q' : ['CAA', 'CAG'],
               'N' : ['AAT', 'AAC'],
               'H' : ['CAT', 'CAC'],
               'E' : ['GAA', 'GAG'],
               'D' : ['GAT', 'GAC'],
               'K' : ['AAA', 'AAG'],
               'R' : ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
               'Stop' : ['TAA', 'TAG', 'TGA'],
				'-' : ['---']}
codon_to_aa = {codon : aa for aa,v in codon_table.items() for codon in v}
