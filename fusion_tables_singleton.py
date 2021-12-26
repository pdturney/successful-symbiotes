#
# Fusion Tables Singleton
#
# Peter Turney, December 15, 2021
#
# Calculate various statistics for a single fusion that has no
# children (a singleton). 
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
# Let N be the maximum number of parts in a symbiote. We need to consider
# values ranging from 0 to N inclusive, so each table will contain N+1 
# rows and columns.
#
max_parts = 5
table_range = max_parts + 1
#
#
# Location of all_seed_storage.bin runs.
#
family_tree_dir = "C:/Users/peter/Peter's Projects" + \
                  "/successful-symbiotes/Experiments"
#
# File for report of fusion table. Use suffix "tsv" for tab-separated
# values, for easy importing of the file into LibreOffice.
#
fusion_table_file = family_tree_dir + "/fusion_tables_singleton.tsv"
fusion_table_handle = open(fusion_table_file, "w")
#
# Make 12 tables for singletons:
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
  # loop through the fusion seeds and look for singletons
  #
  for fusion_seed_ID in fusion_seed_identifiers:
    #
    g.show("run " + str(i+1) + " ... processing seed " + str(fusion_seed_ID) + " ...")
    g.update()
    #
    # here we are only interested in singeltons, so skip over seeds
    # with children
    #
    children = map_parent_to_children[fusion_seed_ID]
    if (len(children) > 0):
      # skip to the next fusion_seed_ID
      continue
    #
    # - number of parts in seed
    #
    fusion_seed = map_ID_to_seed[fusion_seed_ID]
    region_map = mfunc.region_map(fusion_seed)
    num_parts = np.amax(region_map)
    if (num_parts > max_parts):
      # skip to the next fusion_seed_ID
      continue
    #
    # - sample size
    #
    sample_size += 1
    #
    # - run number
    #
    run_number = i + 1
    #
    # - singleton ID
    #
    singleton_ID = fusion_seed_ID
    singleton_seed = fusion_seed
    #
    # - node mutualism stats
    #
    #                [red, blue,   orange,     green]
    colour_weights = [1.0, 0.0, (2.0 / 3.0), (1.0 / 3.0)]
    num_steps = 1000
    [singleton_insider_count, singleton_outsider_count, \
      singleton_insider_growth, singleton_outsider_growth] = \
      mfunc.seed_mutualism_stats(g, singleton_seed, num_steps, colour_weights)
    #
    # - node interaction stats
    #
    [singleton_ensemble_count, singleton_soloist_count, \
      singleton_ensemble_growth, singleton_soloist_growth] = \
      mfunc.seed_interaction_stats(g, singleton_seed, num_steps)
    #
    # - node management stats
    #
    [singleton_manager_count, singleton_worker_count, \
      singleton_manager_growth, singleton_worker_growth] = \
      mfunc.seed_management_stats(g, singleton_seed, num_steps)
    #
    # update the tables
    #
    man_cnt[singleton_manager_count, singleton_worker_count] += 1
    man_sum[singleton_manager_count, singleton_worker_count] += \
      singleton_manager_growth + singleton_worker_growth
    #
    mut_cnt[singleton_outsider_count, singleton_insider_count] += 1
    mut_sum[singleton_outsider_count, singleton_insider_count] += \
      singleton_outsider_growth + singleton_insider_growth
    #
    int_cnt[singleton_soloist_count, singleton_ensemble_count] += 1
    int_sum[singleton_soloist_count, singleton_ensemble_count] += \
      singleton_soloist_growth + singleton_ensemble_growth
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
fusion_table_handle.write("\n\nFusion Tables Singleton\n\n" + \
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