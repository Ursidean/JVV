'''Fix random seed'''
'''
Fix the random number seed in the Metronamica project file to generate multiple
solutions
'''

#Modules
import numpy as np
import xml.etree.ElementTree as ET

#Function
def set_rand(project_path,rseed):

    source=open(project_path)
    tree=ET.parse(source)                                                 #Parse georpoject file as xml type data structure for modification
    root=tree.getroot()
    
    spline_fixed=root[2][1][3][0][0][10][0][0]
    spline_value=root[2][1][3][0][0][10][0][1]
    
    fixed=1

    spline_fixed.text=str(fixed)
    spline_value.text=str(rseed)
    tree.write(project_path)