#
# Fusion Descendants
#
# Peter Turney, November 15, 2021
#
# Read all of the seeds (all_seed_storage.bin) and search for seeds
# that are instances of new fusions (seed.birth_type = "fusion"). Then
# for each new fusion, search for all of its descendants, until the 
# most recent descendant is found. The seeds are numbered by their
# birth order, so the most recent descendant is the descendant with
# the largest number.
#
import golly as g
import model_classes as mclass
import model_functions as mfunc
import model_parameters as mparam
import numpy as np
import copy
import time
import pickle
import os
import re
import sys
#
# Location of all_seed_storage.bin runs.
#
family_tree_dir = "C:/Users/peter/Peter's Projects" + \
                  "/successful-symbiotes/Experiments"
#
# Number of "run" subdirectories.
#
num_runs = 40
#
# Loop through the runs and extract the seeds.
#
for i in range(num_runs):
  #
  g.show("run " + str(i+1) + " ...")
  g.update()
  #
  # read the data from run i and store it in seed_list
  #
  # "Experiments/run1", "Experiments/run2", ...
  seed_file = family_tree_dir + "/run" + str(i+1) + "/all_seed_storage.bin"
  # open the pickle file -- "ab+" opens a file for 
  # both appending and reading in binary mode
  seed_handle = open(seed_file, "ab+")
  # start at the beginning of the file
  seed_handle.seek(0)
  # initialize the list of seeds
  seed_list = []
  # step through all of the pickles in the file, adding them
  # to the list
  while True:
    try:
      seed = pickle.load(seed_handle)
      seed_list.append(seed)
    except (EOFError, pickle.UnpicklingError):
      break
  # close the file
  seed_handle.close()
  #
  g.show("run " + str(i+1) + " ... mapping parents to children ...")
  g.update()
  #
  # each seed specifies its parents (except for the initial
  # random seeds, which have no parents), but we want to know
  # each seed's children, not its parents, so we will step
  # through the list of seeds and construct a mapping from
  # the given seed's ID number to a list of the ID numbers of
  # its children -- that is, in map_parent_to_children, the children
  # lists are listed in the order of their parents in seed_list
  #
  avoid_list = ["random", "fusion", "fission"] # these aren't really "children"
  #
  map_parent_to_children = mfunc.map_parent_to_child(seed_list, avoid_list)
  #
  # make a list of types of seeds, listed in the order given by seed_list
  #
  map_seed_list_to_type = mfunc.find_types(seed_list)
  #
  # look for fusion seeds to use as roots of our trees
  # -- we are only interested in trees that have fusions as their roots
  #
  fusion_seed_identifiers = mfunc.fusion_seed_IDs(seed_list)
  #
  # sample size for randomly sampling full node-to-leaf paths in the tree
  # -- sample_size = 10 ** sample_exponent
  #
  sample_exponent = 6
  exponent_format = "{:." + str(sample_exponent) + "f}"
  #
  # loop through the fusion seeds and create a family tree file for each seed
  #
  for fusion_seed_ID in fusion_seed_identifiers:
    #
    g.show("run " + str(i+1) + " ... processing seed " + str(fusion_seed_ID) + " ...")
    g.update()
    #
    descendants_properties = mfunc.longest_paths(fusion_seed_ID, \
      map_parent_to_children, map_seed_list_to_type, sample_exponent)
    #
    # write out all the full fusion paths to the report file
    #
    report_file = family_tree_dir + "/run" + str(i + 1) + "/fusion" + \
      str(fusion_seed_ID) + "-descendants-sample10power" + str(sample_exponent) + ".txt"
    report_handle = open(report_file, "w")
    #
    # descendants_properties is a sorted list of lists of the following form:
    #
    # [
    #   [node X ID, node X type, node X prob, node X avg depth, node X num children],
    #   [node Y ID, node Y type, node Y prob, node Y avg depth, node Y num children],
    #   ...
    # ]
    for [node_ID, node_type, node_prob, node_avg_depth, node_num_children] \
      in descendants_properties:
      report_handle.write("node: " + str(node_ID) + ",  type: " + str(node_type) + \
        ",  prob: " + exponent_format.format(node_prob) + ",  depth: " + \
        str(node_avg_depth) + ", children: " + str(node_num_children) + "\n")
    #
    report_handle.close()
  #
  g.show("run " + str(i+1) + " done.")
  g.update()
  #
#
#