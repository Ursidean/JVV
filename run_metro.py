'''Run Metronamica'''

#Modules
import os
import subprocess

#Function
def run_metro(project_file,log_file,working_dir,geo_cmd):
    #Move to working directory
    os.chdir(working_dir)
    #First, reset the project file to the first year
    p1=subprocess.Popen([geo_cmd, "--Reset", project_file], shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    for line in p1.stdout.readlines():
        print line
    retval=p1.wait()
    #Second, run the model, generate an output map for analysis
    p2=subprocess.Popen([geo_cmd, "--Run", "--LogSettings", log_file, project_file,], shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    for line in p2.stdout.readlines():
        print line,
    retval=p2.wait()
    
    