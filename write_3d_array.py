''' Three dimensional array printing '''
#Modules
import numpy as np

#Function
def write_3d_array(outfile,data):
    outfile=open(outfile,'w')
    outfile.write('# Array shape: {0}\n'.format(data.shape))
    for data_slice in data:
        np.savetxt(outfile, data_slice, fmt='%-7.2f')
        outfile.write('# New slice\n')
    
    outfile.close()

