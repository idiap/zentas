# Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
# Written by James Newling <jnewling@idiap.ch>
# zentas is a k-medoids library written in C++ and Python. This file is part of zentas.
# zentas is free software: you can redistribute it and/or modify it under the terms of
# the GNU General Public License version 3 as published by the Free Software Foundation.
# zentas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE. See the GNU General Public License for more details. You should have received
# a copy of the GNU General Public License along with zentas. If not, see
# <http://www.gnu.org/licenses/>.

import hardpaths
reload(hardpaths)

import numpy as np
import numpy.random as npr
import hardpaths
reload(hardpaths)

fasta_data_dir = hardpaths.fasta_data_dir
import os

def make_snippet_file(chr_38_fn = os.path.join(fasta_data_dir, "Homo_sapiens.GRCh38.dna.chromosome.10.fa"), n_snips = 10**2 ,snip_range = [20, 30]):
  """
  Read a long string from *chr_38_fn*, and extract *n_snips* snippets of lengths in *snip_range*, and write to *chr_38_fn*.snips.fa
  """
  filly = open(chr_38_fn, "r")
  lines = filly.readlines()
  lines = [x[0:-1] for x in lines[170:]]
  filly.close()
  
  oneline = ''.join(lines)
  
  n_nucleotides = len(oneline)
  print "n nucleotides : ", n_nucleotides
  
  
  snippet_fn = chr_38_fn + ".snips.fa"
  filly = open(snippet_fn, "w")
  
  
  for i in range(n_snips):
    start = npr.randint(0, n_nucleotides - snip_range[-1] - 2)
    end = start + npr.randint(snip_range[0], snip_range[1])
    snippet = oneline[start:end]
    filly.write(">snippet_%d\n%s\n\n"%(i,snippet))
  
  filly.close()
  
