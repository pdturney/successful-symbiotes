#
# Fusion Tables Prolific
#
# Peter Turney, December 15, 2021
#
# Calculate various statistics for a fusion that has children
# (a family of fusions). From each family, select the fusion seed
# that has the most children (a prolific seed).
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
# Number of "run" subdirectories.
#
num_runs = 40
#
# Sample size for randomly sampling full node-to-leaf paths in the tree
# (sample_size = 10 ** sample_exponent).
#
sample_exponent = 6
#
# Let N be the maximum number of parts in a symbiote. We need to consider
# values ranging from 0 to N inclusive, so each table will contain N+1 
# rows and columns.
#
max_parts = 5
table_range = max_parts + 1
#
# Location of all_seed_storage.bin runs.
#
family_tree_dir = "C:/Users/peter/Peter's Projects" + \
                  "/successful-symbiotes/Experiments"
#
# File for report of fusion table. Use suffix "tsv" for tab-separated
# values, for easy importing of the file into LibreOffice.
#
fusion_table_file = family_tree_dir + "/fusion_tables_prolific.tsv"
fusion_table_handle = open(fusion_table_file, "w")
#
# Make 12 tables for prolific nodes:
#
#  (1) Management Count: counts of managers and workers
#  (2) Management Sum: sums of growths of managers and workers
#  (3) Management Average: averages of growths of managers and workers
#  (4) Management Percent: percents of managers and workers
#
#  (5) Mutualism Count: counts of insiders and outsiders
#  (6) Mutualism Sum: sums of growths of insiders and outsiders
#  (7) Mutualism Average: averages of insiders and outsiders
#  (8) Mutualism Percent: percents of insiders and outsiders
#
#  (9) Interaction Count: counts of soloists and ensembles
# (10) Interaction Sum: sums of growths of soloists and ensembles
# (11) Interaction Average: averages of soloists and ensembles
# (12) Interaction Percent: percents of soloists and ensembles
#
man_cnt = np.zeros((table_range, table_range), dtype=np.int)
man_sum = np.zeros((table_range, table_range), dtype=np.int)
man_avg = np.zeros((table_range, table_range), dtype=np.double)
man_pct = np.zeros((table_range, table_range), dtype=np.double)
#
mut_cnt = np.zeros((table_range, table_range), dtype=np.int)
mut_sum = np.zeros((table_range, table_range), dtype=np.int)
mut_avg = np.zeros((table_range, table_range), dtype=np.double)
mut_pct = np.zeros((table_range, table_range), dtype=np.double)
#
int_cnt = np.zeros((table_range, table_range), dtype=np.int)
int_sum = np.zeros((table_range, table_range), dtype=np.int)
int_avg = np.zeros((table_range, table_range), dtype=np.double)
int_pct = np.zeros((table_range, table_range), dtype=np.double)
#
# Loop through the runs and extract the seeds.
#
sample_size = 0
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
  # make a mapping from seed ID numbers to actual seeds
  #
  map_ID_to_seed = {} # dictionary
  for seed in seed_list:
    map_ID_to_seed[seed.unique_ID_num] = seed
  #
  # look for fusion seeds to use as roots of our trees
  # -- we are only interested in trees that have fusions as their roots
  #
  fusion_seed_identifiers = mfunc.fusion_seed_IDs(seed_list)
  #
  # loop through the fusion seeds and create a family tree file for each seed
  #
  for fusion_seed_ID in fusion_seed_identifiers:
    #
    g.show("run " + str(i+1) + " ... processing seed " + str(fusion_seed_ID) + " ...")
    g.update()
    #
    # here we are only interested in families, so skip over seeds
    # with no children
    #
    children = map_parent_to_children[fusion_seed_ID]
    num_children = len(children)
    if (num_children == 0):
      # skip to the next fusion_seed_ID
      continue
    #
    # - number of parts in seed
    #
    # we assume here that all seeds in the same family of descendants
    # must have the same number of parts, so we only need to count the
    # number of parts in the seed that forms the root of the family tree
    #
    fusion_seed = map_ID_to_seed[fusion_seed_ID]
    region_map = mfunc.region_map(fusion_seed)
    num_parts = np.amax(region_map)
    if (num_parts > max_parts):
      # skip to the next fusion_seed_ID
      continue
    #
    # descendants_properties is a sorted list of lists of the following form:
    #
    # [
    #   [node X ID, node X type, node X prob, node X avg depth, node X num children],
    #   [node Y ID, node Y type, node Y prob, node Y avg depth, node Y num children],
    #   ...
    # ]
    #
    descendants_properties = mfunc.longest_paths(fusion_seed_ID, \
      map_parent_to_children, map_seed_list_to_type, sample_exponent)
    #
    # - sample size
    #
    sample_size += 1
    #
    # - run_number
    #
    run_number = i + 1
    #
    # - root_ID (for current fusion family tree)
    #
    root_ID = fusion_seed_ID
    #
    # - last_ID (for current fusion family tree)
    #
    last_ID = descendants_properties[-1][0] # [last row][first column]
    #
    # - prolific_node_ID (ID of seed with most children)
    # - prolific_node_count (count of seed with most children)
    #
    # initialize prolific_node: prolific_node will eventually be
    # updated to have maximum number of children
    prolific_node_ID = fusion_seed_ID
    prolific_node_count = num_children
    # loop through the list of nodes in the tree of descendants
    num_descendants = len(descendants_properties)
    for j in range(num_descendants):
      # relevant properties for node j
      parent_ID = descendants_properties[j][0]
      children_count = descendants_properties[j][4]
      # update prolific_node?
      if (children_count > prolific_node_count):
        prolific_node_ID = parent_ID
        prolific_node_count = children_count
    #
    # - number of nodes in tree (number of descendants)
    #
    num_nodes = len(descendants_properties)
    #
    # - prolific mutualism = [prolific_insider_count, prolific_outsider_count,
    #                         prolific_insider_growth, prolific_outsider_growth]
    #
    #                [red, blue,   orange,     green]
    colour_weights = [1.0, 0.0, (2.0 / 3.0), (1.0 / 3.0)]
    num_steps = 1000
    prolific_node = map_ID_to_seed[prolific_node_ID]
    [prolific_insider_count, prolific_outsider_count, \
      prolific_insider_growth, prolific_outsider_growth] = \
      mfunc.seed_mutualism_stats(g, prolific_node, num_steps, colour_weights)
    #
    # - prolific interaction = [prolific_ensemble_count, prolific_soloist_count,
    #                           prolific_ensemble_growth, prolific_soloist_growth]
    #
    [prolific_ensemble_count, prolific_soloist_count, \
      prolific_ensemble_growth, prolific_soloist_growth] = \
      mfunc.seed_interaction_stats(g, prolific_node, num_steps)
    #
    # - prolific management = [prolific_manager_count, prolific_worker_count,
    #                          prolific_manager_growth, prolific_worker_growth]
    #
    [prolific_manager_count, prolific_worker_count,
      prolific_manager_growth, prolific_worker_growth] = \
      mfunc.seed_management_stats(g, prolific_node, num_steps)
    #
    # update the tables
    #
    man_cnt[prolific_manager_count, prolific_worker_count] += 1
    man_sum[prolific_manager_count, prolific_worker_count] += \
      prolific_manager_growth + prolific_worker_growth
    #
    mut_cnt[prolific_outsider_count, prolific_insider_count] += 1
    mut_sum[prolific_outsider_count, prolific_insider_count] += \
      prolific_outsider_growth + prolific_insider_growth
    #
    int_cnt[prolific_soloist_count, prolific_ensemble_count] += 1
    int_sum[prolific_soloist_count, prolific_ensemble_count] += \
      prolific_soloist_growth + prolific_ensemble_growth
    #
    #
  #
  g.show("run " + str(i+1) + " done.")
  g.update()
  #
#
# calculate the average growths from the sums and the counts
#
for row in range(table_range):
  for col in range(table_range):
    if (man_cnt[row, col] > 0):
      man_avg[row, col] = man_sum[row, col] / man_cnt[row, col]
    else:
      man_avg[row, col] = 0.0
    #
    if (mut_cnt[row, col] > 0):
      mut_avg[row, col] = mut_sum[row, col] / mut_cnt[row, col]
    else:
      mut_avg[row, col] = 0.0
    #
    if (int_cnt[row, col] > 0):
      int_avg[row, col] = int_sum[row, col] / int_cnt[row, col]
    else:
      int_avg[row, col] = 0.0
#
# calculate the percentages from the counts
#
for row in range(table_range):
  for col in range(table_range):
    man_pct[row, col] = man_cnt[row, col] / sample_size
    mut_pct[row, col] = mut_cnt[row, col] / sample_size
    int_pct[row, col] = int_cnt[row, col] / sample_size
#
# write out the tables
#
# - header for tables
#
fusion_table_handle.write("\n\nFusion Tables Prolific\n\n" + \
  "Sample size = " + str(sample_size) + "\n\n")
#
#  (1) Management Count: counts of managers and workers
#
table_matrix = man_cnt
table_handle = fusion_table_handle
table_title = "Management Count: counts of managers and workers"
row_label = "managers"
col_label = "workers"
#
mfunc.write_fusion_tables(table_matrix, table_handle, table_title, \
  table_range, row_label, col_label)
#
#  (2) Management Sum: sums of growths of managers and workers
#
table_matrix = man_sum
table_handle = fusion_table_handle
table_title = "Management Sum: sums of growths of managers and workers"
row_label = "managers"
col_label = "workers"
#
mfunc.write_fusion_tables(table_matrix, table_handle, table_title, \
  table_range, row_label, col_label)
#
#  (3) Management Average: averages of growths of managers and workers
#
table_matrix = man_avg
table_handle = fusion_table_handle
table_title = "Management Average: averages of growths of managers and workers"
row_label = "managers"
col_label = "workers"
#
mfunc.write_fusion_tables(table_matrix, table_handle, table_title, \
  table_range, row_label, col_label)
#
#  (4) Management Percent: percents of managers and workers
#
table_matrix = man_pct
table_handle = fusion_table_handle
table_title = "Management Percent: percents of managers and workers"
row_label = "managers"
col_label = "workers"
#
mfunc.write_fusion_tables(table_matrix, table_handle, table_title, \
  table_range, row_label, col_label)
#
#  (5) Mutualism Count: counts of insiders and outsiders
#
table_matrix = mut_cnt
table_handle = fusion_table_handle
table_title = "Mutualism Count: counts of insiders and outsiders"
row_label = "outsiders"
col_label = "insiders"
#
mfunc.write_fusion_tables(table_matrix, table_handle, table_title, \
  table_range, row_label, col_label)
#
#  (6) Mutualism Sum: sums of growths of insiders and outsiders
#
table_matrix = mut_sum
table_handle = fusion_table_handle
table_title = "Mutualism Sum: sums of growths of insiders and outsiders"
row_label = "outsiders"
col_label = "insiders"
#
mfunc.write_fusion_tables(table_matrix, table_handle, table_title, \
  table_range, row_label, col_label)
#
#  (7) Mutualism Average: averages of insiders and outsiders
#
table_matrix = mut_avg
table_handle = fusion_table_handle
table_title = "Mutualism Average: averages of insiders and outsiders"
row_label = "outsiders"
col_label = "insiders"
#
mfunc.write_fusion_tables(table_matrix, table_handle, table_title, \
  table_range, row_label, col_label)
#
#  (8) Mutualism Percent: percents of insiders and outsiders
#
table_matrix = mut_pct
table_handle = fusion_table_handle
table_title = "Mutualism Percent: percents of insiders and outsiders"
row_label = "outsiders"
col_label = "insiders"
#
mfunc.write_fusion_tables(table_matrix, table_handle, table_title, \
  table_range, row_label, col_label)
#
#  (9) Interaction Count: counts of soloists and ensembles
#
table_matrix = int_cnt
table_handle = fusion_table_handle
table_title = "Interaction Count: counts of soloists and ensembles"
row_label = "soloists"
col_label = "ensembles"
#
mfunc.write_fusion_tables(table_matrix, table_handle, table_title, \
  table_range, row_label, col_label)
#
# (10) Interaction Sum: sums of growths of soloists and ensembles
#
table_matrix = int_sum
table_handle = fusion_table_handle
table_title = "Interaction Sum: sums of growths of soloists and ensembles"
row_label = "soloists"
col_label = "ensembles"
#
mfunc.write_fusion_tables(table_matrix, table_handle, table_title, \
  table_range, row_label, col_label)
#
# (11) Interaction Average: averages of soloists and ensembles
#
table_matrix = int_avg
table_handle = fusion_table_handle
table_title = "Interaction Average: averages of soloists and ensembles"
row_label = "soloists"
col_label = "ensembles"
#
mfunc.write_fusion_tables(table_matrix, table_handle, table_title, \
  table_range, row_label, col_label)
#
# (12) Interaction Percent: percents of soloists and ensembles
#
table_matrix = int_pct
table_handle = fusion_table_handle
table_title = "Interaction Percent: percents of soloists and ensembles"
row_label = "soloists"
col_label = "ensembles"
#
mfunc.write_fusion_tables(table_matrix, table_handle, table_title, \
  table_range, row_label, col_label)
#
# close report file
#
fusion_table_handle.close()
#
#