#!/usr/bin/python -u

#########################################################################
#########################################################################
##                                                                     ##
##             PYCAL: 2-dimensional heat equation solver               ##
##   for arbitary geometries (structured and unstructured grids).      ##
##                                                                     ##
#########################################################################
##                                                                     ##
##  Mesh needs to be generated with gmsh (http://geuz.org/gmsh/)       ##
##                                                                     ##
#########################################################################
#########################################################################


import argparse
import numpy as np
import scipy as sp
import math
import sys
import os.path
import time
import fvm_elem as msh

        
class input_data:
    def __init__(self):
        self.BCT = []
        self.BC = []
        self.h_coeffs = []
        self.T0 = 298.0
        self.t_total = 1000.0
        self.dt = 100.0
        ## The default properties correspond to air at 25C and 100KPa (SI units).
        self.k_cond = 0.0257
        self.rho = 1.1839
        self.cp = 1.005
        
    def read_input(self, ff):
        for line in ff:
            if line.split('=')[0] == "BCT " or line.split('=')[0] == "BCT":
                for i in range(len(line.split('=')[1].split(','))):
                    self.BCT.append(int(line.split('=')[1].split(',')[i]))
            if line.split('=')[0] == "BC " or line.split('=')[0] == "BC":
                for i in range(len(line.split('=')[1].split(','))):
                    self.BC.append(float(line.split('=')[1].split(',')[i]))
            if line.split('=')[0] == "h_coeffs " or line.split('=')[0] == "h_coeffs":
                for i in range(len(line.split('=')[1].split(','))):
                    self.h_coeffs.append(float(line.split('=')[1].split(',')[i]))
            if line.split('=')[0] == "T0 " or line.split('=')[0] == "T0":
                self.T0 = float(line.split('=')[1])
            if line.split('=')[0] == "t_total " or line.split('=')[0] == "t_total":
                self.t_total = float(line.split('=')[1])
            if line.split('=')[0] == "k_cond " or line.split('=')[0] == "k_cond":
                self.k_cond = float(line.split('=')[1])
            if line.split('=')[0] == "rho " or line.split('=')[0] == "rho":
                self.rho = float(line.split('=')[1])
            if line.split('=')[0] == "cp " or line.split('=')[0] == "cp":
                self.cp = float(line.split('=')[1])
            if line.split('=')[0] == "dt " or line.split('=')[0] == "dt":
                self.dt = float(line.split('=')[1])
    


#########################################################################
##                                                                     ##
##  Main program starts here.                                          ##
##                                                                     ##
#########################################################################


start_time = time.time()


parser = argparse.ArgumentParser(description='2-D Heat equation solver for structured and unstructured meshes.')
parser.add_argument('-i','--input', help='Input file name (plain text)',required=True) 
parser.add_argument('-m','--mesh', help='Mesh input file name (.msh)',required=True)
parser.add_argument('--restart', help='Restart Calculation.', action='store_true',required=False)


args = parser.parse_args()
meshFile = args.mesh
inputFile = args.input
outputFile=meshFile.split('inp')[0]+"out"+meshFile.split('inp')[1].split('.msh')[0]+".out"
restartSimulation = args.restart


stopFile=meshFile.split('.msh')[0]+".stop"


if os.path.exists(stopFile):
    print "Please delete "+stopFile+" before proceeding."
    quit()

if restartSimulation:
    if not os.path.exists(outputFile):
        print "There is no previous output file."
        quit()


##  Input data is read from the input file and stored in input_data object

f=open(inputFile,'r')

data_in = input_data()
data_in.read_input(f)

f.close()

# Time divisions
M = int(data_in.t_total/data_in.dt)

# Thermal diffusivity (m2/s)
alpha = data_in.k_cond/(data_in.rho*data_in.cp)


##  Mesh data is read from the mesh file and stored in Node and Element objects

f=open(meshFile,'r')
meshData=[]
for line in f:
    meshData.append(line)


for i in range(len(meshData)):
    if meshData[i]=='$Nodes\n':
        N_nodes = int(meshData[i+1])
        begin_nodes = i+2
        break

for i in range(len(meshData)):
    if meshData[i]=='$Elements\n':
        N_elements = int(meshData[i+1])
        begin_elements = i+2
        break

f.close()

       
meshNodes = []
for i in range(begin_nodes,begin_nodes+N_nodes):
    aux_n = msh.Nodes(int(meshData[i].split(' ')[0]), float(meshData[i].split(' ')[1]), float(meshData[i].split(' ')[2]))
    meshNodes.append(aux_n)

meshElements = []
non_elements = 1
for i in range(begin_elements,begin_elements+N_elements):
    aux_id = int(meshData[i].split(' ')[0])
    aux_eletype = int(meshData[i].split(' ')[1])
    aux_flag = int(meshData[i].split(' ')[2])
    aux_bdid = int(meshData[i].split(' ')[3])
    aux_nodes = []
    if aux_eletype == 1:
        aux_n_nodes = 2
        for j in range(5,5+aux_n_nodes):
            aux_nodes.append(meshNodes[int(meshData[i].split(' ')[j])-1])
    if aux_eletype == 2:
        aux_n_nodes = 3
        for j in range(5,5+aux_n_nodes):
            aux_nodes.append(meshNodes[int(meshData[i].split(' ')[j])-1])
    if aux_eletype == 3:
        aux_n_nodes = 4
        for j in range(5,5+aux_n_nodes):
            aux_nodes.append(meshNodes[int(meshData[i].split(' ')[j])-1])
    if aux_eletype != 15:
        aux_elements = msh.Elements(aux_id, aux_eletype, aux_flag, aux_bdid, aux_n_nodes, aux_nodes)
        meshElements.append(aux_elements)
    else:
        non_elements = non_elements+1


print "Number of elements: "+str(len(meshElements))
print "Number of timesteps: "+str(M)
print "Total simulation time: "+str(data_in.t_total/60)+" min"

 

##  Element calculations are performed.

for i in range(len(meshElements)):
    if meshElements[i].elementType == 1:
        for j in range(len(data_in.BCT)):
            if meshElements[i].boundaryId == 200+j:
                if data_in.BCT[j] == 1:
                    meshElements[i].boundaryFlux = data_in.h_coeffs[j]*data_in.BC[j]
                elif data_in.BCT[j] == 2:
                    meshElements[i].boundaryFlux = data_in.BC[j]

total_area = 0.0
for i in range(len(meshElements)):
    for j in range(i,len(meshElements)):
        meshElements[i].find_neighbours(meshElements[j])
        if len(meshElements[i].neighbours) == 4:
            break
    meshElements[i].calc_normals()
    meshElements[i].calc_edgeMidpoint_centroid()
    meshElements[i].calc_area()
    total_area = total_area + meshElements[i].area



for i in range(len(meshElements)):
    for j in range(len(meshElements[i].neighbours)):
#        meshElements[i].calc_lambda(meshElements[meshElements[i].neighbours[j]-non_elements],j)
        meshElements[i].calc_dx(meshElements[meshElements[i].neighbours[j]-non_elements])
        meshElements[i].grad_dir.append(meshElements[i].centroid.diff_vec(meshElements[meshElements[i].neighbours[j]-non_elements].centroid))
        meshElements[i].grad_dir[j].normalise()
        if meshElements[i].elementType == 2 or meshElements[i].elementType == 3:
            if meshElements[i].norms[j].dot_prod(meshElements[i].MidpointToCentroid[j]) > 0:
              meshElements[i].norms[j].invert_vec()


               
# Assembly of Matrices

N = 0
for i in range(len(meshElements)):
    N = N+1

matrixA = np.zeros(shape=(N,N))
matrixB = np.zeros(shape=(N,N))
boundary_flux = np.zeros(shape=(N,1))

Temps = np.zeros(shape=(N,M+1))


for i in range(len(meshElements)):
    Temps[i][0] = data_in.T0
    if meshElements[i].elementType == 2 or meshElements[i].elementType == 3:
        F = []
        for j in range(len(meshElements[i].neighbours)):
            a = meshElements[i].norms[j].dot_prod(meshElements[i].grad_dir[j])
            F.append(-(meshElements[i].edges[j].norm)/(meshElements[i].area*meshElements[i].ds[j]))
        Flux_total = 0.0
        for j in range(len(meshElements[i].neighbours)):
            Flux_total = Flux_total + F[j]
        matrixA[i][i] = 1.0 - (alpha*data_in.dt*Flux_total)
        for j in range(len(meshElements[i].neighbours)):
            matrixA[i][meshElements[i].neighbours[j]-1] = alpha*data_in.dt*F[j]
        matrixB[i][i] = 1.0
    if meshElements[i].elementType == 1:
        for j in range(len(meshElements[meshElements[i].neighbours[0]-1].neighbours)):
            if meshElements[meshElements[meshElements[i].neighbours[0]-1].neighbours[j]-1].elementType == 1:
                PointToPoint = meshElements[meshElements[i].neighbours[0]-1].MidpointToCentroid[j].norm
        for j in range(len(data_in.BCT)):
            if meshElements[i].boundaryId == 200+j:
                if data_in.BCT[j] == 1:
                    matrixA[i][i] = -1.0
                    matrixA[i][meshElements[i].neighbours[0]-1] = 1.0
                elif data_in.BCT[j] == 2:
                    matrixA[i][i] = 1.0
        for j in range(len(data_in.BCT)):
            if meshElements[i].boundaryId == 200+j:
                if data_in.BCT[j] == 1:
                    matrixB[i][i] = data_in.h_coeffs[j]*PointToPoint/data_in.k_cond
                    boundary_flux[i] = -meshElements[i].boundaryFlux*PointToPoint/data_in.k_cond
                elif data_in.BCT[j] == 2:
                    boundary_flux[i] = meshElements[i].boundaryFlux

A = np.asmatrix(matrixA)
B = np.asmatrix(matrixB)
Tboundary = np.asmatrix(boundary_flux)


# The initial conditions for restarting a simulation are set.

start = 0
replaceStart = 0 
aux = 0
out_old = []

if restartSimulation:
    g = open(outputFile,'r')
    for line in g:
        if line == "$ElementData\n":
            aux = 1
        if aux == 1:
           out_old.append(line) 
    if int(out_old[8]) != len(meshElements):
        print "The number of elements does not match!"
        quit()

    n = len(meshElements)+4
    M_old = 2*int(out_old[-n])
    if M_old >= M:
        print "Previous simulation shorter than new simulation."
        quit()
    else:
        start = M_old
    ii = 0
    print start
    replaceStart = int(out_old[-n])
    print replaceStart
    for j in range((replaceStart)*(len(meshElements)+10),(replaceStart+1)*(len(meshElements)+10)-1):
        if ii > 8:
            Temps[ii-9][start] = float(out_old[j].split(' ')[1])
        ii = ii +1
    g.close()
    g = open(outputFile,'a')
    replaceStart = replaceStart + 1
else:
    g=open(outputFile,'w')
    for i in range(len(meshData)):
        g.write(meshData[i])


# Solution iterations

for i in range(start,M):
    if os.path.exists(stopFile):
        M = i
        break
    print "Time step "+str(i)
    Tn = np.asmatrix(Temps[:,[i]])
    # Tnpo = (A.I*B)*Tn + A.I*Tboundary
    Tnpo = np.linalg.solve(matrixA, np.dot(matrixB,Tn) + boundary_flux) 
    Temps[:,[i+1]] = Tnpo
    



# Output is written to the output file

ii = replaceStart
for i in xrange(start,M,1):
    g.write("$ElementData\n")
    g.write("1\n")
    g.write("\"Temperature [K]\"\n")
    g.write("1\n")
    g.write(str(i*data_in.dt)+"\n")
    g.write("3\n")
    g.write(str(ii)+"\n")
    ii = ii + 1
    g.write("1\n")
    g.write(str(N)+"\n")
    for j in range(N):
        g.write(str(meshElements[j].e_id)+" "+str(Temps[j][i])+"\n")

    g.write("$EndElementData\n")

g.close()

print "It took "+str(time.time()-start_time)+" seconds to run the simulation."

