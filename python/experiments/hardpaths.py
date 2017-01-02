# Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
# Written by James Newling <jnewling@idiap.ch>
# zentas is a k-medoids library written in C++ and Python. This file is part of zentas. zentas is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. zentas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with zentas. If not, see <http://www.gnu.org/licenses/>.

import os
import socket
import sys

if (socket.gethostname() == 'idbean'):
  raise RuntimeError("You need to set some paths if running on idbean")

else:   #(socket.gethostname() == 'goudurix12'):
  username = "j" + "n" + "e" + "w" + "l" + "i" + "n" + "g"
  institute = "i" + "d" + "i" + "a" + "p"
  
  homedir = os.path.join("/" , institute, "home", username)
  tempdir = os.path.join("/" , institute, "temp", username)
  userdir = os.path.join("/" , institute, "user", username)

fasta_data_dir = os.path.join(tempdir, "data/fastadata")
infipath = os.path.join(tempdir,  "data/infinitemnist/infimnist")
infiexec = os.path.join(infipath, "infimnist")
zentas_dir = os.path.join(homedir, "libraries/zentasactive")
zentas_data_dir = os.path.join(zentas_dir, "data")
zentas_output_dir = os.path.join(userdir, "zentasoutput")
clarans_vik_results_base = os.path.join(zentas_output_dir, "results/clarans_vik_experiments/")
english_words_dir = os.path.join(userdir, "english-words")
fasta_data_dir = os.path.join(tempdir,"data/fastadata/")
initialisation_experiments_runtimes_dir = os.path.join(zentas_output_dir, "initialisation_experiments_runtimes/")
initialisation_result_pickles = os.path.join(zentas_output_dir, "initialisation_result_pickles")
elapsed_time_kmeanspp = os.path.join(initialisation_result_pickles, "elapsed_time_kmeanspp")
clarans_paper_dir = os.path.join(homedir, "colftex/clarans/")
zentas_lib_dir = "../build/python"
trainandtestdir = os.path.join(tempdir, "data/densedata/trainandtest")
joensuudata = os.path.join(userdir, "joensuudata")
