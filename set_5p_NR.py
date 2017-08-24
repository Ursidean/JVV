'''Adjust neighbourhood rules'''
'''
Used to set neighbourhood rules in Metronamica project file. Neighbourhood rules
are assumed to be structured as follows (5 point stance):
|
o
|
|   o
|  
|     o  
|
|           o    
|
|                   o
|                   
|         
---------------------------------o-----

Where the circles are:
The influence at 0 (inertia or conversion)
The influence at 1 
The influence at 2
The location where influence is 0
'''

#Modules
import numpy as np
import xml.etree.ElementTree as ET

#Function
def set_5p_rule(project_path,fu_elem,lu_elem,y0,y1,y2,y3,y4):
    #y0 is the influence at a distance of 0
    #y1 is the influence at a distance of 1
    #y2 is the influence at a distance of 2
    #y3 is the influence at a distance of 3
    #y4 is the influence at a distance of 4
    source=open(project_path)
    tree=ET.parse(source)                                                 #Parse georpoject file as xml type data structure for modification
    root=tree.getroot()
    #Access specified neighbourhood rule
    spline=root[2][1][3][0][0][1][0][1][0][0][fu_elem][0][lu_elem][0]
    #Store array of values
    arrange=np.array([[0,y0],[1,y1],[2,y2],[3,y3],[4,y4],[5,0]])
    #Find length of arrange
    size4=len(arrange)
    #Remove current neighbourhood rule values
    for point in spline.findall('point'):
        spline.remove(point)
    #Dictionary to store inputs for neighbourhood rules
    inputs={}                                                             
    for i in range(0,size4):
        key='new_value'+str(i)
        inputs[key]={'y':str(arrange[i,1]),'x':str(arrange[i,0])}
    #Insert points, commit to case study    
    for i in range(0,size4):                                              
        key='new_value'+str(i)
        child=ET.SubElement(spline, 'point', attrib=inputs[key]) 
    
    tree.write(project_path)
    source.close()