"""
Validation of neighbourhood rules obtained from the implementation of Van Vliet
et al. (2013) method across both calibration and validation time periods.
"""

# Allows for multi-dimensional array handling.
import numpy as np
# Read input maps as 2d arrays.
from read_map import read_map
# Find the size of different square neighbourhoods.
from considered_distances import considered_distances
# Calculate teh enrichment factor.
from enrichment_factor import ef
# Set the random seed.
from set_rand import set_rand
# Set the neighoburhood rule parameter values.
from set_5p_NR import set_5p_rule
# Run Metronamica through the command line.
from run_metro import run_metro
# Convert numbers to log base 10
from math import log10
# Read and write from csv file
import csv
# Function to calculate Kappa and Kappa Simulation.
from kappa import kappa
from kappa import ksim
# Function to calculate the average area weighted clumpiness error.
from area_weighted_clu import area_weighted_clu_error

# Specify the base path of hte directory.
base_path = "C:\\Users\\charl\\OneDrive\\Documents\\JVV\\"
# Specify the case study application.
case_study = "Madrid"
# Set the paths to directories containing data and output.
data_path = base_path + "EU_test_data\\"
output_path = base_path + "EU_test_output\\" + case_study
# Specify the map paths.
map1_path = data_path + case_study + "\\" + case_study.lower() + "_1990.asc"
map2_path = data_path + case_study + "\\" + case_study.lower() + "_2000.asc"
map3_path = data_path + case_study + "\\" + case_study.lower() + "_2006.asc"
mask_path = data_path + case_study + "\\" + case_study.lower() + "_mask.asc"
# Specify the working directory for the calibration data.
working_directory_cal = ("C:\\Geonamica\\Metronamica\\" +
                         case_study)
# Specify the working directory for the validation data.
working_directory_val = working_directory_cal + "_2000"
# Specify the project file for the calibraiton data.
project_file_cal = working_directory_cal + "\\" + case_study + ".geoproj"
# Specify the project file for the validation data.
project_file_val = working_directory_val + "\\" + case_study + "_2000.geoproj"
# Specify the path to the simulated output map for the calibration period.
smap_path_cal = (working_directory_cal + "\\Log\\Land_use\\" 
                 "Land use map_2000-Jan-01 00_00_00.rst")
# Specify the path to the simulated output map for the validation period.
smap_path_val = (working_directory_val + "\\Log\\Land_use\\" 
                 "Land use map_2006-Jan-01 00_00_00.rst")
# Specify the log file for the year 2000.
log_file_2000 = base_path + "LogSettings.xml"
log_file_2006 = base_path + "LogSettings_2006.xml"
# Specify the command line version of Geonamica.
geo_cmd = "C:\\Program Files (x86)\\Geonamica\\Metronamica\\GeonamicaCmd.exe"
# Set the land-use class names.
luc_names = ["Natural areas", "Arable land", "Permanent crops", "Pastures",
             "Agricultural areas", "Residential", "Industry & commerce",
             "Recreation areas", "Forest", "Road & rail", "Seaports", 
             "Airports", "Mine & dump sites", "Fresh water", "Marine water"]
# Set the land-use class parameters: number of land-use classes, passive, 
# feature and active.
luc = len(luc_names)
pas = 1
fea = 6
act = luc - (pas + fea)
# Read in the map for the data at time slice 1.
map_1990 = read_map(map1_path)
# Read in the map for the data at time slice 2.
map_2000 = read_map(map2_path)
# Read in the map for the data at time slice 3.
map_2006 = read_map(map3_path)
# Read in the masking map
mask = read_map(mask_path)
# Analyse the input maps for evaluation purposes.
map_dimensions = np.shape(map_1990)
rows = map_dimensions[0]
cols = map_dimensions[1]
# Specify the maximum neighoburhood size distance considered.
max_distance = 5
# Determine the distances that will be analysed.
temp = considered_distances(max_distance)
# Store the list of considered distances as a variable.
cd = temp[0]
# Store the total number of distances considered.
cdl = temp[1]
# Determine the maximum neighbourhood size (unit) from considered distances.
N_all = [1, 8, 12, 16, 32, 28, 40, 40, 20]
N = []
for c in range(0, max_distance):
    N.append(N_all[c])
# Count the presence of each land-use class for the calibration and validation
# data maps.
luc_count_2000 = [0]*luc
luc_count_2006 = [0]*luc
for i in range(0, rows):
    for j in range(0, cols):
        if mask[i, j] > 0:
            luc_count_2000[map_2000[i, j]] = luc_count_2000[map_2000[i, j]] + 1
            luc_count_2006[map_2006[i, j]] = luc_count_2006[map_2006[i, j]] + 1
    
    
# Determine the enrichment factor values for the data for the calibration time 
# period.
data_ef_cal = ef(luc, max_distance, cdl, cd, N, map_1990, map_2000, mask, rows, 
                 cols)
# Determine the enrichment factor values for the data for the validation time 
# period.
data_ef_val = ef(luc, max_distance, cdl, cd, N, map_2000, map_2006, mask, rows, 
                 cols)
# Specify the input rules file.
input_rule_file = output_path + "\\Rules\\JVV_rules.csv"
# Read the input rule files.
rules = {}
for i in range(0, luc):
    for j in range(0, act):
        key = "from " + luc_names[i] + " to " + luc_names[j + pas]
        rules[key] = [0]*5
# Read inputs from a csv file.
with open (input_rule_file, 'rb') as f:
    readCSV = csv.reader(f)
    # Skip the header line.
    next(f)
    for row in readCSV:
        i = row[0]
        j = row[1]
        key = "from " + i + " to " + j
        rules[key][0] = row[2]
        rules[key][1] = row[3]
        rules[key][2] = row[4]
        rules[key][3] = row[5]
        rules[key][4] = row[6]
# Input the rule set into the models.
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
        set_5p_rule(project_file_cal, fu_elem, lu_elem, y0, y1, y2, y3, y4)
        set_5p_rule(project_file_val, fu_elem, lu_elem, y0, y1, y2, y3, y4)
        
# Evaluate the simulated output for the calibrated and validated output for the 
# set number of evaluations.
evaluations = 5
base_seed = 1000
# Initialise a set of lists to store values.
fuzzy_kappa_cal = [0]*evaluations
fuzzy_kappa_val = [0]*evaluations
fks_cal = [0]*evaluations
fks_val = [0]*evaluations
clu_cal = [0]*evaluations
clu_val = [0]*evaluations
max_dev_cal = [0]*evaluations
max_dev_val = [0]*evaluations
# Set the constant value for the calculation of deviation.
c = 0.25
# Begin the iterative testing.
for x in range(0, evaluations):
    # Set the random seed in the models.
    rseed = base_seed + x
    set_rand(project_file_cal, rseed)
    set_rand(project_file_val, rseed)
    # Run the model for the calibration period and read the input map.
    run_metro(project_file_cal, log_file_2000, working_directory_cal, geo_cmd)
    smap_2000 = read_map(smap_path_cal)
    # Calculate the fuzzy kappa value.
    fuzzy_kappa_cal[x] = kappa(map_2000, smap_2000, mask)
    # Calculate the fuzzy kappa simulation value.
    fks_cal[x] = ksim(map_1990, map_2000, smap_2000, mask)
    # Calculate the average area weighted clumpiness error
    clu_cal[x] =  area_weighted_clu_error(map_2000, smap_2000, mask, luc, pas, 
                                          act, luc_count_2000)
    # Determine the maximum deviation.
    smap_2000_ef = ef(luc, max_distance, cdl, cd, N, map_1990, smap_2000, mask,
                      rows, cols)
    dev = np.zeros(shape=(max_distance, luc, luc))
    for d in range(0, max_distance):
        for p in range(0, luc):
            for q in range(0, luc):
                dev[d, p, q] = log10((smap_2000_ef[d, p, q] + c) /
                                     (data_ef_cal[d, p, q] + c))
    max_dev_temp = 0
    # Identify the largest meaningful deviation.
    for d in range(0, max_distance):
        for p in range(0, act):
            for q in range(0, luc):
                # Skip conversions for feature classes.
                if q > (act + pas - 1) and d == 0:
                    pass
                else:
                    if abs(dev[d, p + pas, q]) > max_dev_temp:
                        max_dev_temp = abs(dev[d, p + pas, q])
    # Assign value to list.
    max_dev_cal[x] = max_dev_temp
    
    # Now perform analysis for the validation period.
    # Run the model for the validation period and read the input map.
    run_metro(project_file_val, log_file_2006, working_directory_val, geo_cmd)
    smap_2006 = read_map(smap_path_val)
    # Calculate the fuzzy kappa value.
    fuzzy_kappa_val[x] = kappa(map_2006, smap_2006, mask)
    # Calculate the fuzzy kappa simulation value.
    fks_val[x] = ksim(map_2000, map_2006, smap_2006, mask)
    # Calculate the average area weighted clumpiness error
    clu_val[x] =  area_weighted_clu_error(map_2006, smap_2006, mask, luc, pas, 
                                          act, luc_count_2006)
    # Determine the maximum deviation.
    smap_2006_ef = ef(luc, max_distance, cdl, cd, N, map_2000, smap_2006, mask,
                      rows, cols)
    dev = np.zeros(shape=(max_distance, luc, luc))
    for d in range(0, max_distance):
        for p in range(0, luc):
            for q in range(0, luc):
                dev[d, p, q] = log10((smap_2006_ef[d, p, q] + c) /
                                     (data_ef_val[d, p, q] + c))
    max_dev_temp = 0                                
    # Identify the largest meaningful deviation.
    for d in range(0, max_distance):
        for p in range(0, act):
            for q in range(0, luc):
                # Skip conversions for feature classes.
                if q > (act + pas - 1) and d == 0:
                    pass
                else:
                    if abs(dev[d, p + pas, q]) > max_dev_temp:
                        max_dev_temp = abs(dev[d, p + pas, q])
    # Assign value to list.
    max_dev_val[x] = max_dev_temp                         

# Write the output to a csv file.
metrics_output_file = output_path + "\\metrics.csv"
store = [0] * 5
with open (metrics_output_file, "wb") as csv_file:
    writer = csv.writer(csv_file)
    # Save and write a header line to the file.
    values = ["Period", "Fuzzy Kappa", "FKS", "CLU error", "Max dev."]
    writer.writerow(values)
    for x in range(0, evaluations):
        store[0] = "Cal"
        store[1] = fuzzy_kappa_cal[x]
        store[2] = fks_cal[x]
        store[3] = clu_cal[x]
        store[4] = max_dev_cal[x]
        writer.writerow(store)
    for x in range(0, evaluations):
        store[0] = "Val"
        store[1] = fuzzy_kappa_val[x]
        store[2] = fks_val[x]
        store[3] = clu_val[x]
        store[4] = max_dev_val[x]
        writer.writerow(store)

# Finished!
