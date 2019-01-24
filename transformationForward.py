'''
author : hashem
05/04/2018
'''

import sys, os
import numpy as np
import math
import random
import json
# Intialise OpenCMISS/iron 
from opencmiss.iron import iron

rotation = True
if (rotation):
    translation = False
else: 
    translation = True

#=================================================================
# Control Panel
#=================================================================
# set the number of elements and the number of nodes for the cylinder 
numberOfDimensions = 3
numberOfGaussXi = 3 
numberOfCircumfrentialElementsPerQuarter = 2
numberOfCircumfrentialElements = 4*numberOfCircumfrentialElementsPerQuarter
numberOfCircumfrentialNodes = numberOfCircumfrentialElements
numberOfLengthElements = 8
numberOfLengthNodes = numberOfLengthElements+1
numberOfWallElements = 1
numberOfWallNodes = numberOfWallElements+1
origin = [0.0,0.0,0.0]
meshOrigin = [0.0,0.0,0.0]


numberOfNodes = numberOfCircumfrentialElements*(numberOfLengthElements+1)*(numberOfWallElements+1)
if (numberOfWallElements == 0): 
    numberOfElements = numberOfCircumfrentialElements*numberOfLengthElements
elif (numberOfWallElements == 1): 
    numberOfElements = numberOfCircumfrentialElements*numberOfLengthElements
else:
    numberOfElements = numberOfCircumfrentialElements*numberOfLengthElements*numberOfWallElements


# the results of the optimisation
#optmin = [21.44, -3.21, -12.82, -80.27, -164.54]
#optmin = [ 0.0, 0.0, 0.0]
#optmin = [  ]

#optmin = [ -6.17 , 6.42 , 12.70 ]
#optmin = [-0.4826, -0.00769, -0.0769 ]
#optmin = [-80.0, 6.6, 12.0 ]
#optmin = [19.30, -6.23, 6.23, -173.11, -144.65, 1.96 ]
#optmin = [-0.288, -0.128,-0.128, -59, -67, 13.70]
#optmin = [20.4203520709, -12.4052172626, -7.29874303596, -186.376191643, -249.214856552, -3.99027091452]
#optmin = [0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
optmin = [ -0.234948271993, -0.0148436611266, 0.0148433232846]
stage = 20
initialLocation = np.zeros((numberOfNodes, 3))
with open("Coords20.txt", "r") as ins:
    arrayOfInitialInputData = []
    for line in ins:
        arrayOfInitialInputData.append(line)
x,y,z = 0.0,0.0,0.0
for i in range (numberOfNodes):
    for j in range (3):
        sample = arrayOfInitialInputData[i*3 + j]
        if (math.fmod(j,3) == 0):
            x = float (sample[0:25])                
        elif (math.fmod(j,3) == 1):
            y = float (sample[0:25])
        elif (math.fmod(j,3) == 2):
            z = float (sample[0:25])
            initialLocation[i,:] = [x,y,z]

finalLocation = np.zeros((numberOfNodes, 3))
with open("Coords0.txt", "r") as ins:
    arrayOfInputData = []
    for line in ins:
        arrayOfInputData.append(line)
x,y,z = 0.0,0.0,0.0
for i in range (numberOfNodes):
    for j in range (3):
        sample = arrayOfInputData[i*3 + j]
        if (math.fmod(j,3) == 0):
            x = float (sample[0:25])                
        elif (math.fmod(j,3) == 1):
            y = float (sample[0:25])
        elif (math.fmod(j,3) == 2):
            z = float (sample[0:25])
            finalLocation[i,:] = [x,y,z]

originFirst = np.zeros((1,3))
# bringing the data of the meshes to the origin  
originFirst[0,:] = initialLocation[139 - 1 , :]
for i in range (numberOfNodes):
    initialLocation [i,:] = initialLocation [i,:] - originFirst[0,:]

originSecond = np.zeros((1,3))
originSecond[0,:] = finalLocation[139 - 1 , :]
for i in range (numberOfNodes):
    finalLocation [i,:] = finalLocation [i,:] - originSecond[0,:]


#print initialLocation
#print originSecond
#sys.exit()

# =========================================================
#
#         Setting the optimizer for the problem 
#          Finding the transformation function 
#
# =========================================================
numberOfXcomponents = 3
numberOfYcomponents = 3
translationComponents = 3
#There are three different sets of optimization needed to work for this problem 
#  the for of the transformation matrix is: 
#
#               | u(0) , u(1) , u(2) , [u(9) ]| 
#          T =  | u(3) , u(4) , u(5) , [u(10)]| 
#               | u(6) , u(7) , u(8) , [u(11)]|                   
#               |   0  ,   0  ,  0   ,   1    |     
#
# defining the objective function for the  optimization Number one 
# probably the rates are symmetric ... 


aa = math.cos(optmin[0])*math.cos(optmin[1]) - math.cos(optmin[1])*math.sin(optmin[2])*math.sin(optmin[0])             #   cos(Psi) * cos(Theta) - cos(Theta) * sin(Phi) * sin(Psi)
ab = math.cos(optmin[0])*math.sin(optmin[2]) + math.cos(optmin[1])*math.cos(optmin[2])*math.sin(optmin[0])             #   cos(Psi) * sin(Phi)   + cos(Theta) * cos(Phi) * sin(Psi)
ac = math.sin(optmin[0])*math.sin(optmin[1])                                                                           #   sin(Psi) * sin(Theta)

ba =-math.sin(optmin[0])*math.cos(optmin[2]) - math.cos(optmin[1])*math.sin(optmin[2])*math.cos(optmin[0])             # - sin(Psi) * cos(Phi)   - cos(Theta) * sin(Phi) * cos(Psi) 
bb =-math.sin(optmin[2])*math.sin(optmin[0]) + math.cos(optmin[1])*math.cos(optmin[2])*math.cos(optmin[0])             # - sin(Phi) * sin(Psi)   - cos(Theta) * cos(Phi) * cos(Psi) 
bc = math.cos(optmin[0])*math.sin(optmin[1])                                                                           #   cos(Psi) * sin(Theta)

ca = math.sin(optmin[1])*math.sin(optmin[2])                                                                           #   sin(Theta) * sin(Phi)
cb =-math.sin(optmin[1])*math.cos(optmin[2])                                                                           # - sin(Theta) * cos(Phi)
cc = math.cos(optmin[1])                                                                                               #   cos (Theta)

T = np.zeros((3, 3))
V = np.zeros((3,1))
T[0,:] =aa,ab,ac
T[1,:] =ba,bb,bc
T[2,:] =ca,cb,cc


det = np.linalg.det(T)
print det 
if (rotation):
    print "pure rotation"
else: 
    V[0,0] = optmin [3]
    V[1,0] = optmin [4]
    V[2,0] = 0.0

print V

vect = np.zeros((3))
newFinalLocation = np.zeros((numberOfNodes, 3))
for i in range (numberOfNodes):
    vect[:] = initialLocation[i,:]
    if (rotation):
        newFinalLocation[i,:] = T.dot(vect)
    else:
        newFinalLocation[i,:] = T.dot(vect)  + V[:,0]

#print newInitialLocation

#=================================================================
# User Numbers , CS, Region, Basis
#=================================================================
(coordinateSystemUserNumber,
    regionUserNumber,
    basisUserNumber,
    generatedMeshUserNumber,
    meshUserNumber,
    decompositionUserNumber,
    geometricFieldUserNumber,
    equationsSetFieldUserNumber,
    dependentFieldUserNumber,
    independentFieldUserNumber,
    dataPointFieldUserNumber,
    materialFieldUserNumber,
    analyticFieldUserNumber,
    dependentDataFieldUserNumber,
    dataPointsUserNumber,
    dataProjectionUserNumber,
    equationsSetUserNumber,
    problemUserNumber) = range(1,19)
numberOfComputationalNodes = iron.ComputationalNumberOfNodesGet()
computationalNodeNumber = iron.ComputationalNodeNumberGet()
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.dimension = 3
coordinateSystem.CreateFinish()
region = iron.Region()
region.CreateStart(regionUserNumber,iron.WorldRegion)
region.label = "FittingRegion"
region.coordinateSystem = coordinateSystem
region.CreateFinish()
basis = iron.Basis()
basis.CreateStart(basisUserNumber)
basis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
basis.numberOfXi = numberOfDimensions
basis.interpolationXi = [iron.BasisInterpolationSpecifications.CUBIC_HERMITE]*3
basis.quadratureNumberOfGaussXi = [numberOfGaussXi]*3
basis.CreateFinish()

#=================================================================
# Mesh & Decomposition
#=================================================================
# creating the number of elements and the mesh origins ... and/or
# Start the creation of a manually generated mesh in the region
numberOfNodes = numberOfCircumfrentialElements*(numberOfLengthElements+1)*(numberOfWallElements+1)
if (numberOfWallElements == 0): 
    numberOfElements = numberOfCircumfrentialElements*numberOfLengthElements
elif (numberOfWallElements == 1): 
    numberOfElements = numberOfCircumfrentialElements*numberOfLengthElements
else:
    numberOfElements = numberOfCircumfrentialElements*numberOfLengthElements*numberOfWallElements

mesh = iron.Mesh()
mesh.CreateStart(meshUserNumber,region,3)
mesh.origin = meshOrigin
mesh.NumberOfComponentsSet(1)
mesh.NumberOfElementsSet(numberOfElements)
# Define nodes for the mesh
nodes = iron.Nodes()
nodes.CreateStart(region,numberOfNodes)
nodes.CreateFinish()
elements = iron.MeshElements()
meshComponentNumber = 1
elements.CreateStart(mesh, meshComponentNumber, basis)
elementNumber = 0
for wallElementIdx in range(1,numberOfWallElements+1):
    for lengthElementIdx in range(1,numberOfLengthElements+1):
        for circumfrentialElementIdx in range(1,numberOfCircumfrentialElements+1):
            elementNumber = elementNumber + 1
            localNode1 = circumfrentialElementIdx + (lengthElementIdx - 1)*numberOfCircumfrentialElements + \
                (wallElementIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes
            if circumfrentialElementIdx == numberOfCircumfrentialElements:
                localNode2 = 1 + (lengthElementIdx-1)*numberOfCircumfrentialNodes + \
                    (wallElementIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes
            else: 
                localNode2 = localNode1 + 1
            localNode3 = localNode1 + numberOfCircumfrentialNodes
            localNode4 = localNode2 + numberOfCircumfrentialNodes
            localNode5 = localNode1 + numberOfCircumfrentialNodes*numberOfLengthNodes
            localNode6 = localNode2 + numberOfCircumfrentialNodes*numberOfLengthNodes
            localNode7 = localNode3 + numberOfCircumfrentialNodes*numberOfLengthNodes
            localNode8 = localNode4 + numberOfCircumfrentialNodes*numberOfLengthNodes
            localNodes = [localNode1,localNode2,localNode3,localNode4,localNode5,localNode6,localNode7,localNode8]
            elements.NodesSet(elementNumber,localNodes)  
elements.CreateFinish()
mesh.CreateFinish() 
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber,mesh)
decomposition.type = iron.DecompositionTypes.CALCULATED
decomposition.numberOfDomains = numberOfComputationalNodes
decomposition.CalculateFacesSet(True)
decomposition.CreateFinish()
#  ===============================================================
#  =============================================================== reading stage 10 
manualNodePoints = np.zeros((numberOfLengthNodes,numberOfCircumfrentialNodes,3,numberOfWallNodes))


for wallNodeIdx in range(1,numberOfWallNodes+1):
    for lengthNodeIdx in range(1,numberOfLengthNodes+1):
        for circumfrentialNodeIdx in range(1,numberOfCircumfrentialNodes+1):
            nodeNumber = circumfrentialNodeIdx + (lengthNodeIdx-1)*numberOfCircumfrentialNodes + (wallNodeIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes 
            manualNodePoints[lengthNodeIdx-1,circumfrentialNodeIdx-1,:,wallNodeIdx-1] = newFinalLocation[nodeNumber-1,:]


#print manualNodePoints
manualPoints = np.zeros((numberOfLengthNodes*numberOfCircumfrentialNodes*numberOfWallNodes,3,8))
referencePoint = np.zeros((3,1))

#referencePoint [0,0] = manualNodePoints[8,2,0,0] 
#referencePoint [1,0] = manualNodePoints[8,2,1,0] 
#referencePoint [2,0] = manualNodePoints[8,2,2,0] 

#for wallNodeIdx in range(1,numberOfWallNodes+1):
#    for lengthNodeIdx in range(1,numberOfLengthNodes+1):
#        for circumfrentialNodeIdx in range(1,numberOfCircumfrentialNodes+1):
#            manualNodePoints[lengthNodeIdx-1,circumfrentialNodeIdx-1,:,wallNodeIdx-1] = manualNodePoints[lengthNodeIdx-1,circumfrentialNodeIdx-1,:,wallNodeIdx-1] - referencePoint [:,0]

# forming the derivatives with the average of ... 
#calculating the derivatives 
difference = np.zeros((numberOfLengthNodes,numberOfCircumfrentialNodes,3,2))
differenceAverage = np.zeros((numberOfLengthNodes,numberOfCircumfrentialNodes,3,2))
circumDeriv = np.zeros((numberOfLengthNodes,numberOfCircumfrentialNodes,3,2))
directDeriv = np.zeros((numberOfLengthNodes,numberOfCircumfrentialNodes,3,2))
lengthDeriv = np.zeros((numberOfLengthNodes,numberOfCircumfrentialNodes,3,2))
#circumferential derivative to be calculated 
for k in range (2):
    for j in range (numberOfLengthNodes):
        for i in range (numberOfCircumfrentialNodes):
            if (i<numberOfCircumfrentialNodes-1):
                for m in range (3):
                    difference[j,i,m,k]=manualNodePoints[j,i+1,m,k]-manualNodePoints[j,i,m,k]
            else:
                for m in range (3):
                    difference[j,i,m,k]=manualNodePoints[j,0,m,k]-manualNodePoints[j,numberOfCircumfrentialNodes-1,m,k]
for k in range (2):
    for j in range (numberOfLengthNodes):
        for i in range (numberOfCircumfrentialNodes):
            if (i<numberOfCircumfrentialNodes-1):
                for m in range (3):
                    differenceAverage[j,i+1,m,k]=(difference[j,i+1,m,k]+difference[j,i,m,k])/2
            else:
                for m in range (3):
                    differenceAverage[j,0,m,k]=(difference[j,0,m,k]+difference[j,numberOfCircumfrentialNodes-1,m,k])/2
for k in range (2):
    for j in range (numberOfLengthNodes):
        for i in range (numberOfCircumfrentialNodes):
            for m in range (3):
                circumDeriv[j,i,m,k]=differenceAverage[j,i,m,k]/math.sqrt(math.pow(differenceAverage[j,i,0,k],2) + math.pow(differenceAverage[j,i,1,k],2) + math.pow(differenceAverage[j,i,2,k],2))
# derivative of the length direction
for k in range (2):
    for i in range (numberOfCircumfrentialNodes):
        for j in range (numberOfLengthNodes):
            if (j<numberOfLengthNodes-1):
                for m in range (3):
                    difference[j,i,m,k]=manualNodePoints[j+1,i,m,k]-manualNodePoints[j,i,m,k]
            else:
                for m in range (3):
                    difference[j,i,m,k]=manualNodePoints[j,i,m,k]-manualNodePoints[j-1,i,m,k]
for k in range (2):
    for i in range (numberOfCircumfrentialNodes):
        for j in range (numberOfLengthNodes):
            if (j == 0):
                for m in range (3): 
                    differenceAverage[j,i,m,k]=difference[j,i,m,k]
            if (j<numberOfLengthNodes-1):
                for m in range (3):
                    differenceAverage[j+1,i,m,k]=(difference[j,i,m,k]+difference[j+1,i,m,k])/2
            else:
                for m in range (3):
                    differenceAverage[j,i,m,k]=difference[j-1,i,m,k]
for k in range (2):
    for j in range (numberOfLengthNodes):
        for i in range (numberOfCircumfrentialNodes):
            for m in range (3):
                lengthDeriv[j,i,m,k]=differenceAverage[j,i,m,k]/math.sqrt(math.pow(differenceAverage[j,i,0,k],2) + math.pow(differenceAverage[j,i,1,k],2) + math.pow(differenceAverage[j,i,2,k],2))
# the derivatives of the wall direction is defined in the below lines ... 
for i in range (numberOfCircumfrentialNodes):
    for j in range (numberOfLengthNodes):
        for m in range (3):
            for k in range (2):
                difference[j,i,m,k] = manualNodePoints[j,i,m,1] - manualNodePoints[j,i,m,0]
for i in range (numberOfCircumfrentialNodes):
    for j in range (numberOfLengthNodes):
        for k in range (2):
            for m in range (3):
                directDeriv[j,i,m,k] = difference[j,i,m,k]/math.sqrt(math.pow(difference[j,i,0,k],2) + math.pow(difference[j,i,1,k],2) + math.pow(difference[j,i,2,k],2))

# Create a field for the geometry
geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber,region)
geometricField.meshDecomposition = decomposition
for dimension in range(3):
    geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,dimension+1,1)
geometricField.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
geometricField.CreateFinish()

# Get nodes
nodes = iron.Nodes()
region.NodesGet(nodes)
numberOfNodes = nodes.numberOfNodes
# Create the geometric field
for wallNodeIdx in range(1,numberOfWallNodes+1):
    for lengthNodeIdx in range(1,numberOfLengthNodes+1):
        for circumfrentialNodeIdx in range(1,numberOfCircumfrentialNodes+1):
            nodeNumber = circumfrentialNodeIdx + (lengthNodeIdx-1)*numberOfCircumfrentialNodes + (wallNodeIdx-1)*numberOfCircumfrentialNodes*numberOfLengthNodes 
            x = manualNodePoints[lengthNodeIdx-1, circumfrentialNodeIdx-1, 0, wallNodeIdx-1]
            y = manualNodePoints[lengthNodeIdx-1, circumfrentialNodeIdx-1, 1, wallNodeIdx-1]
            z = manualNodePoints[lengthNodeIdx-1, circumfrentialNodeIdx-1, 2, wallNodeIdx-1]
            xtangent = circumDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 0, wallNodeIdx-1]
            ytangent = circumDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 1, wallNodeIdx-1]
            ztangent = circumDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 2, wallNodeIdx-1]
            xnormal = directDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 0, wallNodeIdx-1]
            ynormal = directDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 1, wallNodeIdx-1]
            znormal = directDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 2, wallNodeIdx-1]
            zxnormal = lengthDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 0, wallNodeIdx-1]
            zynormal = lengthDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 1, wallNodeIdx-1]
            zznormal = lengthDeriv[lengthNodeIdx-1, circumfrentialNodeIdx-1, 2, wallNodeIdx-1]
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,1,nodeNumber,1,x)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,1,nodeNumber,2,y)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,1,nodeNumber,3,z)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,1,xtangent)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,2,ytangent)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S1,nodeNumber,3,ztangent)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,1,zxnormal)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,2,zynormal)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S2,nodeNumber,3,zznormal)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,1,xnormal)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,2,ynormal)
            geometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                        1,iron.GlobalDerivativeConstants.GLOBAL_DERIV_S3,nodeNumber,3,znormal)
    # Update the geometric field
geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
    # Export undeformed mesh geometry
######    print("Writing undeformed geometry")
fields = iron.Fields()
fields.CreateRegion(region)
fields.NodesExport("newMesh"+str(stage)+"-8x8","FORTRAN")
fields.ElementsExport("newMesh"+str(stage)+"-8x8","FORTRAN")
fields.Finalise()

iron.Finalise()

