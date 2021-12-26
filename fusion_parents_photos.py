#
# Fusion Parents Photos
#
# Peter Turney, November 14, 2021
#
# Read all of the seeds (all_seed_storage.bin) and search for seeds
# that are instances of new fusions (seed.birth_type = "fusion"). Then,
# for each new fusion, search for its two parents (each fusion must
# have exactly two parents). Take a photo of each fusion and its two
# parents.
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
seed_dir = "C:/Users/peter/Peter's Projects" + \
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
  seed_file = seed_dir + "/run" + str(i+1) + "/all_seed_storage.bin"
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
  # For each fusion seed, we will take one photo of the initial
  # state of the fusion seed and then one photo of the initial
  # state of each of the two parent seeds. Thus we'll have three
  # photos for each fusion seed.
  #
  photo_directory = seed_dir + "/run" + str(i+1)
  # extract leaf directory from photo_directory (so we know where it came from, 
  # in case it gets moved)
  leaf_dir = os.path.basename(os.path.normpath(photo_directory))
  # allow time for Golly image to stabilize before entering loop below
  time.sleep(2)
  # pause between images, in seconds
  pause = 0.1
  #
  # step through the seeds, looking for fusion seeds, and note
  # the parents of the fusion seeds
  #
  for seed in seed_list:
    if (seed.birth_type == "fusion"):
      # get the fusion seed identifier
      seed_ID = seed.unique_ID_num
      # show some info in Golly
      g.show("processing seed " + str(seed_ID) + "...")
      g.update()
      # get some info about the parents
      parent_A_ID = seed.parent_A_ID_num
      parent_B_ID = seed.parent_B_ID_num
      parent_A = seed_list[parent_A_ID]
      parent_B = seed_list[parent_B_ID]
      parent_A_type = parent_A.birth_type
      parent_B_type = parent_B.birth_type
      #
      # photo 1: the fusion seed
      #
      steps = 0 # initial state
      file_path = photo_directory + "/" + leaf_dir + "-fusion-seed" + \
                  str(seed_ID) + ".png"
      rule_name = "Management"
      description = "fusion seed ID " + str(seed_ID)
      mfunc.snap_photo(g, file_path, rule_name, seed, \
                       steps, description, pause)
      #
      # photo 2: parent A
      #
      steps = 0 # initial state
      file_path = photo_directory + "/" + leaf_dir + "-fusion-seed" + \
                  str(seed_ID) + "-parent-A" + str(parent_A_ID) + ".png"
      rule_name = "Management"
      description = "fusion seed ID " + str(parent_A_ID)
      mfunc.snap_photo(g, file_path, rule_name, parent_A, \
                       steps, description, pause)
      #
      # photo 3: parent B
      #
      steps = 0 # initial state
      file_path = photo_directory + "/" + leaf_dir + "-fusion-seed" + \
                  str(seed_ID) + "-parent-B" + str(parent_B_ID) + ".png"
      rule_name = "Management"
      description = "fusion seed ID " + str(parent_B_ID)
      mfunc.snap_photo(g, file_path, rule_name, parent_B, \
                       steps, description, pause)
      #
      # write out the maps of the regions as text files
      # for inspection and verification
      #
      file_path = photo_directory + "/" + leaf_dir + "-fusion-seed" + \
                  str(seed_ID) + "-region-maps.txt"
      file_handle = open(file_path, "w")
      file_handle.write("\n\nRegion map for fusion seed and parents\n\n")
      # process the fusion seed and its two parents
      for seed in [seed, parent_A, parent_B]:
        file_handle.write("Seed " + str(seed.unique_ID_num) + "\n\n")
        seed_map = mfunc.region_map(seed)
        num_rows = seed_map.shape[0]
        num_cols = seed_map.shape[1]
        # rotate matrix for consistency with Golly
        for y in range(num_cols): 
          for x in range(num_rows):
            file_handle.write("{:4d}".format(seed_map[x][y]))
          file_handle.write("\n")
        file_handle.write("\n\n")
      file_handle.close()
      #
    #
  #
#
#