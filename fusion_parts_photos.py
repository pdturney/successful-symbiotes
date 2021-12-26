#
# Fusion Parts Photos
#
# Peter Turney, November 14, 2021
#
# Read fusion pickle files (fusion_storage.bin) and
# make photos of all possible pairwise colourings, where
# one part is red and all other parts are blue.
#
# It is expected that this code will be executed inside Golly
# and that the Golly window will be maximized, so that the photos
# will have maximum resolution.
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
# Parameters
#
num_steps = 1000 # number of steps to run Management Rule
num_files = 40   # number of "run" folders: "run1", "run2", ...
#
# Location of fusion_storage.bin files -- the input pickles.
#
fusion_dir = "C:/Users/peter/Peter's Projects" + \
             "/successful-symbiotes/Experiments"
#
# Loop through the fusion files and extract the seeds.
#
for i in range(num_files):
  #
  # assume the folders have the form "run1", "run2", ...
  fusion_file = fusion_dir + "/run" + str(i+1) + "/fusion_storage.bin"
  # open the fusion pickle file -- "ab+" opens a file for 
  # both appending and reading in binary mode
  fusion_handle = open(fusion_file, "ab+")
  # start at the beginning of the file
  fusion_handle.seek(0)
  # initialize the list of seeds
  fusion_list = []
  # step through all of the pickles in the file, adding them
  # to the list
  while True:
    try:
      item = pickle.load(fusion_handle)
      fusion_list.append(item)
    except (EOFError, pickle.UnpicklingError):
      break
  #
  fusion_handle.close()
  #
  # The list fusion_list is a repeating sequence of four items:
  #
  # [s2, s3, s4, n, ..., s2, s3, s4, n]
  #
  # - s2 is part of s4 (after rotation)
  # - s3 is part of s4 (after rotation)
  # - s4 is the fusion of s2 and s3
  # - s4 is the n-th child born
  #
  # For each [s2, s3, s4, n] tuple, if s4 has p parts, then we
  # will compare each part to the union of all other parts. The
  # focal part will be red and all other parts will be blue.
  #
  # For each of the p parts, we will take one photo of the initial
  # state and one photo of the final state. Thus there will be
  # a total of p photos.
  #
  photo_directory = fusion_dir + "/run" + str(i+1)
  # extract leaf directory from photo_directory (so we know where it came from, 
  # in case it gets moved)
  leaf_dir = os.path.basename(os.path.normpath(photo_directory))
  # allow time for Golly image to stabilize before entering loop below
  time.sleep(2)
  # pause between images, in seconds
  pause = 0.1
  #
  # read four items at a time -- we only care about s4 and n, but
  # the fusion pickle file is read by several other scripts, and some
  # of these scripts use s2 and s3 -- s4 is the whole seed and n is the 
  # birth number of s4
  #
  for (s2, s3, s4, n) in zip(*[iter(fusion_list)] * 4):
    #
    # make a matrix that is a map of the regions in s4; the map will
    # have the same size matrix as s4 and each cell in the map will 
    # contain a number that uniquely defines the region it belongs to,
    # as determined by the purple borders in s4
    #
    s4_map = mfunc.region_map(s4)
    #
    # adjust the birth number: in previous versions, we didn't count
    # the births of the initial random population; now we do
    #
    revised_birth_num = n + mparam.pop_size
    #
    # write the map out as a text file for inspection and verification
    #
    file_path = photo_directory + "/" + leaf_dir + "-birth" + \
                str(revised_birth_num) + "-region-map.txt"
    num_rows = s4_map.shape[0]
    num_cols = s4_map.shape[1]
    file_handle = open(file_path, "w")
    # rotate matrix for consistency with Golly
    for y in range(num_cols): 
      for x in range(num_rows):
        file_handle.write("{:4d}".format(s4_map[x][y]))
      file_handle.write("\n")
    file_handle.close()
    #
    # find out how many different regions there are in the map
    # and then generate the p 1-vs-(p-1) colourings
    #
    num_regions = np.amax(s4_map)
    #
    for target_region in range(1, num_regions + 1):
      # make a copy of s4 -- we will create a new colouring
      # for the copy
      s4_colouring = copy.deepcopy(s4)
      # - if a cell in the region map s4_map has the value target_region
      #   and the corresponding cell in s4_colouring is 1 or 2 (red or blue),
      #   then the corresponding cell in s4_colouring will be set
      #   to red (state 1)
      # - if a cell in the region map s4_map does not have the value target_region
      #   and the corresponding cell in s4_colouring is 1 or 2 (red or blue),
      #   then the corresponding cell in s4_colouring will be set
      #   to blue (state 2)
      # - if a cell in the region map s4_map has the value -1,
      #   then the corresponding cell in s4_colouring should be
      #   purple (state 5 -- the colour of the borders); if it is not
      #   purple, then signal an error
      for x in range(num_rows): 
        for y in range(num_cols):
          # state 1 -- red -- target_region
          if ((s4_map[x][y] == target_region) and \
             ((s4_colouring.cells[x][y] == 1) or \
             (s4_colouring.cells[x][y] == 2))):
            s4_colouring.cells[x][y] = 1
          # state 2 -- blue -- not target_region
          elif ((s4_map[x][y] != target_region) and \
               ((s4_colouring.cells[x][y] == 1) or \
               (s4_colouring.cells[x][y] == 2))):
            s4_colouring.cells[x][y] = 2
          # state 5 -- purple -- the border between regions
          elif (s4_map[x][y] == -1):
            assert s4_colouring.cells[x][y] == 5
      #
      # now take two photos of s4_colouring, one for the intial state
      # and another for the final state
      #
      # photo 1: a photo of s4 in its initial state 
      #
      steps = 0 # initial state
      file_path = photo_directory + "/" + leaf_dir + "-birth" + \
                  str(revised_birth_num) + "-target" + str(target_region) + \
                  "-steps" + str(steps) + ".png"
      rule_name = "Management"
      description = "child number " + str(revised_birth_num) + \
                    "initial state, Management"
      mfunc.snap_photo(g, file_path, rule_name, s4_colouring, \
                       steps, description, pause)
      #
      # photo 2: a photo s4 in its final state 
      #
      steps = num_steps # final state
      file_path = photo_directory + "/" + leaf_dir + "-birth" + \
                  str(revised_birth_num) + "-target" + str(target_region) + \
                  "-steps" + str(steps) + ".png"
      rule_name = "Management"
      description = "child number " + str(revised_birth_num) + \
                    "final state, Management"
      mfunc.snap_photo(g, file_path, rule_name, s4_colouring, \
                       steps, description, pause)
      #
      #
    #
    #
  #
  #
#
#