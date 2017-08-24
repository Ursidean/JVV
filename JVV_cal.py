"""
Implementation of Van Vliet et al. (2013) method for the automatic
calibration of neighbourhood rules.
"""

import numpy as np
from read_map import read_map
from considered_distances import considered_distances
from enrichment_factor import ef
from set_rand import set_rand
from set_5p_NR import set_5p_rule
from run_metro import run_metro
from math import log10
import csv

# Specify the base path of the directory
base_path = "C:\\Users\\charl\\OneDrive\\Documents\\JVV\\"
# Select an example case study application. Specify the name below:
case_study = "Budapest"
## Set the paths to the directories and relevant data
data_path = base_path + "Case_study_data\\"
output_path = base_path + "Output\\" + case_study
map1_path = data_path + case_study + "\\" + case_study.lower() + "_1990.asc"
map2_path = data_path + case_study + "\\" + case_study.lower() + "_2000.asc"
mask_path = data_path + case_study + "\\" + case_study.lower() + "_mask.asc"
# Specify the working directory (location of project file) and project file.
working_directory = ("C:\\Geonamica\\Metronamica\\" + case_study)
project_file = (working_directory + "\\" + case_study + ".geoproj")
# Specify the paths to the simulated output maps.
smap_path = (working_directory + "\\Log\\Land_use\\"
                                 "Land use map_2000-Jan-01 00_00_00.rst")
log_file = base_path + "LogSettings.xml"
# Specify the command line version of Geonamica.
geo_cmd = "C:\\Program Files (x86)\\Geonamica\\Metronamica\\GeonamicaCmd.exe"
# Set the land-use class names.
luc_names = ["Natural areas", "Arable land", "Permanent crops", "Pastures",
             "Agricultural areas", "Residential", "Industry & commerce",
             "Recreation areas", "Forest", "Road & rail", "Seaports",
             "Airports", "Mine & dump sites", "Fresh water", "Marine water"]
# Set the land-use class parameters: number of land-use classes, passive,
# feature, and active.
luc = len(luc_names)
pas = 1
fea = 6
act = luc - (pas + fea)
# Read in the map for the data at time slice 1.
omap = read_map(map1_path)
# Read in the map for the data at time slice 2.
amap = read_map(map2_path)
# Read in the masking map.
mask = read_map(mask_path)
# Analyse the input maps for evaluation purposes
map_dimensions = np.shape(omap)
rows = map_dimensions[0]
cols = map_dimensions[1]
# Specify the maximum neighbourhood size distance considered
max_distance = 5
# Determine the distances that will be analysed, use module: considered_distances.
temp = considered_distances(max_distance)
# Store the list of considered distances as a variable.
cd = temp[0]
# Store the total number of distances considered
cdl = temp[1]
# Determine the maximum neighbourhood size (unit) from considered distances
N_all = [1, 8, 12, 16, 32, 28, 40, 40, 20]
N = []
for c in range(0, max_distance):
    N.append(N_all[c])


# Determine the enrichment factor values for the data.
data_ef = ef(luc, max_distance, cdl, cd, N, omap, amap, mask, rows, cols)

# Generate the null rule set
rules = {}
for i in range(0, act):
    for j in range(0, luc):
        key = "from " + luc_names[j] + " to " + luc_names[i + pas]
        if i + pas == j:
            rules[key] = [100, 0, 0, 0, 0]
        else:
            rules[key] = [1, 0, 0, 0, 0]
# Input the null rule set into the model.
for i in range(0, luc):
    for j in range(0, act):
        key = "from " + luc_names[i] + " to " + luc_names[j + pas]
        fu_elem = j
        lu_elem = i
        y0 = rules[key][0]
        y1 = rules[key][1]
        y2 = rules[key][2]
        y3 = rules[key][3]
        y4 = rules[key][4]
        set_5p_rule(project_file, fu_elem, lu_elem, y0, y1, y2, y3, y4)
# Fix the random seed.
rseed = 1000
set_rand(project_file, rseed)
# Set the fixed parameters for iterative testing.
counter = 0
total_iterations = 5000
c = 0.25
# Initialise a matrix to track rule adjustment
rule_tracker=np.zeros(shape=(total_iterations,5))
# Begin the iterative testing.
while counter < total_iterations:
    # Run model to generate a simulated output map.
    run_metro(project_file, log_file, working_directory, geo_cmd)
    # Read the simulation map
    smap = read_map(smap_path)
    # Determine the enrichment factor values for the simualted output map.
    smap_ef = ef(luc, max_distance, cdl, cd, N, omap, smap, mask, rows, cols)
    # Generate the deviation matrix.
    dev = np.zeros(shape=(max_distance, luc, luc))
    for d in range(0, max_distance):
        for p in range(0, luc):
            for q in range(0, luc):
                dev[d, p, q] = log10((smap_ef[d, p, q] + c) /
                                     (data_ef[d, p, q] + c))
    # Initialise a metric to track the total deviation.
    total_dev = 0
    # Initialise variables to find the index of the largest deviation.
    index = [0, 0, 0]
    max_value = 0
    sign = 1
    # Identify the largest deviation.
    for d in range(0, max_distance):
        for p in range(0, act):
            for q in range(0, luc):
                if q > (act + pas - 1) and d == 0:
                    # Skip conversion points for feature classes.
                    # These have no impact on the model.
                    pass
                else:
                    total_dev = total_dev + abs(dev[d, p + pas, q])
                    if abs(dev[d, p + pas, q]) > max_value:
                        index = [d, p + pas, q]
                        max_value = abs(dev[d, p + pas, q])
                        if dev[d, p + pas, q] > 0:
                            sign = 1
                        else:
                            sign = -1
    # Set the adjustment distance.
    adj_d = index[0]
    # Set the adjustment value.
    adj_value = 10**(-1*adj_d)
    # Adjust the corresponding neighbourhood rule in the dictionary.
    rule_key = "from " + luc_names[index[2]] + " to " + luc_names[index[1]]
    if sign == 1:
        rules[rule_key][adj_d] = rules[rule_key][adj_d] - adj_value
    else:
        rules[rule_key][adj_d] = rules[rule_key][adj_d] + adj_value
    # Adjust rule_key value in the Metronamica project file.
    lu_elem = index[2]
    fu_elem = index[1] - pas
    y0 = rules[rule_key][0]
    y1 = rules[rule_key][1]
    y2 = rules[rule_key][2]
    y3 = rules[rule_key][3]
    y4 = rules[rule_key][4]
    set_5p_rule(project_file, fu_elem, lu_elem, y0, y1, y2, y3, y4)
    # Update tracking information.
    rule_tracker[counter, 0] = lu_elem
    rule_tracker[counter, 1] = fu_elem
    rule_tracker[counter, 2] = adj_d
    rule_tracker[counter, 3] = rules[rule_key][adj_d]
    rule_tracker[counter, 4] = total_dev

    # Add one to iteration counter, prevent infinite loop.
    counter = counter + 1
    # User feedback
    print "Iterations completed: " + str(counter)

# Write the output rules.
output_rules_file = output_path + "\\Rules\\JVV_rules.csv"
store = [0]*7
with open(output_rules_file, "wb") as csv_file:
    writer = csv.writer(csv_file)
    # Write a header line to the file
    values = ["from", "to", "y0", "y1", "y2", "y3", "y4"]
    writer.writerow(values)
    # Now write the neighbourhood rules in the form from ... to ...
    for i in range(0, luc):
        for j in range(0, act):
            key = "from " + luc_names[i] + " to " + luc_names[j + pas]
            store[0] = luc_names[i]
            store[1] = luc_names[j + pas]
            store[2] = rules[key][0]
            store[3] = rules[key][1]
            store[4] = rules[key][2]
            store[5] = rules[key][3]
            store[6] = rules[key][4]
            writer.writerow(store)
# Write the output log.
log_file = output_path + "\\JVV_output.csv"
store = [0]*5
with open(log_file, "wb") as csv_file:
    writer = csv.writer(csv_file)
    # Write a header line to the file
    values = ["from", "to", "distance", "new value", "total deviation"]
    writer.writerow(values)
    # Now write the output
    for i in range(0, total_iterations):
        store[0] = luc_names[int(rule_tracker[i, 0])]
        store[1] = luc_names[int(rule_tracker[i, 1] + pas)]
        store[2] = rule_tracker[i, 2]
        store[3] = rule_tracker[i, 3]
        store[4] = rule_tracker[i, 4]
        writer.writerow(store)

#Completed!
