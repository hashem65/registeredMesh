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



def findObjective(optmin,**kwargs):                                           # optmin[0] = Psi  ,   optmin[1] = Theta   ,  optmin[2] = Phi 
    cmissobject=kwargs['cmiss']
    numberOfWallNodes= cmissobject['numberOfWallNodes']
    numberOfCircumfrentialNodes= cmissobject['numberOfCircumfrentialNodes']
    g = []
    func = 0 
    aa = math.cos(optmin[0])*math.cos(optmin[1]) - math.cos(optmin[1])*math.sin(optmin[2])*math.sin(optmin[0])             #   cos(Psi) * cos(Theta) - cos(Theta) * sin(Phi) * sin(Psi)
    ab = math.cos(optmin[0])*math.sin(optmin[2]) + math.cos(optmin[1])*math.cos(optmin[2])*math.sin(optmin[0])             #   cos(Psi) * sin(Phi)   + cos(Theta) * cos(Phi) * sin(Psi)
    ac = math.sin(optmin[0])*math.sin(optmin[1])                                                                           #   sin(Psi) * sin(Theta)

    ba =-math.sin(optmin[0])*math.cos(optmin[2]) - math.cos(optmin[1])*math.sin(optmin[2])*math.cos(optmin[0])             # - sin(Psi) * cos(Phi)   - cos(Theta) * sin(Phi) * cos(Psi) 
    bb =-math.sin(optmin[2])*math.sin(optmin[0]) + math.cos(optmin[1])*math.cos(optmin[2])*math.cos(optmin[0])             # - sin(Phi) * sin(Psi)   - cos(Theta) * cos(Phi) * cos(Psi) 
    bc = math.cos(optmin[0])*math.sin(optmin[1])                                                                           #   cos(Psi) * sin(Theta)

    ca = math.sin(optmin[1])*math.sin(optmin[2])                                                                           #   sin(Theta) * sin(Phi)
    cb =-math.sin(optmin[1])*math.cos(optmin[2])                                                                           # - sin(Theta) * cos(Phi)
    cc = math.cos(optmin[1])                                                                                               #   cos (Theta)

    T = np.array([aa,ab,ac,ba,bb,bc,ca,cb,cc]).reshape((3,3))
     
    det = np.linalg.det(T)
#    if (det != 1):
         
#        return 10e20,g,1

#    aa = 2*(optmin[0]*optmin[0]+optmin[3]*optmin[3]) - 1
#    ab = 2*(optmin[0]*optmin[1]-optmin[2]*optmin[3]) 
#    ac = 2*(optmin[0]*optmin[2]+optmin[1]*optmin[3]) 

#    ba = 2*(optmin[0]*optmin[1]+optmin[2]*optmin[3]) 
#    bb = 2*(optmin[1]*optmin[1]+optmin[3]*optmin[3]) - 1
#    bc = 2*(optmin[1]*optmin[2]-optmin[0]*optmin[3]) 

#    ca = 2*(optmin[0]*optmin[2]-optmin[1]*optmin[3]) 
#    cb = 2*(optmin[1]*optmin[2]+optmin[0]*optmin[3]) 
#    cc = 2*(optmin[2]*optmin[2]+optmin[3]*optmin[3]) - 1

    xlist = [83 ] #,91, 131 ]# ,99,107,115,123,131]
    ylist = [73,77]
    xylist = [74,76,82,84] 
    zlist = [143]
    originlist = [139]

    if (rotation):
        for i in xlist:	
            func+=2*(finalLocation[i,0]-(aa*initialLocation[i,0]+ ab*initialLocation[i,1]+ ac*initialLocation[i,2]))**2
        for i in ylist:	
            func+=(finalLocation[i,1]-(ba*initialLocation[i,0]+ bb*initialLocation[i,1]+ bc*initialLocation[i,2]))**2
        for i in zlist:
            func+=(finalLocation[i,2]-(ca*initialLocation[i,0]+ cb*initialLocation[i,1]+ cc*initialLocation[i,2]))**2
       # for i in xylist:
        #    func+=((finalLocation[i,0]-(aa*initialLocation[i,0]+ ab*initialLocation[i,1]+ ac*initialLocation[i,2]))**2)/2
        #    func+=((finalLocation[i,1]-(ba*initialLocation[i,0]+ bb*initialLocation[i,1]+ bc*initialLocation[i,2]))**2)/2   
    else:
        for i in xlist:	
            func+=2*(finalLocation[i,0]-(aa*initialLocation[i,0]+ ab*initialLocation[i,1]+ ac*initialLocation[i,2]+ optmin[3]))**2
        for i in ylist:	
            func+=(finalLocation[i,1]-(ba*initialLocation[i,0]+ bb*initialLocation[i,1]+ bc*initialLocation[i,2]+ optmin[4]))**2
        for i in zlist:
            func+=(finalLocation[i,2]-(ca*initialLocation[i,0]+ cb*initialLocation[i,1]+ cc*initialLocation[i,2]+ optmin[5]))**2
#        for i in xylist:
#            func+=((finalLocation[i,0]-(aa*initialLocation[i,0]+ ab*initialLocation[i,1]+ ac*initialLocation[i,2]+ optmin[3]))**2)/2
#            func+=((finalLocation[i,1]-(ba*initialLocation[i,0]+ bb*initialLocation[i,1]+ bc*initialLocation[i,2]+ optmin[4]))**2)/2   
        for i in originlist:
            func+= ((finalLocation[i,0]-(aa*initialLocation[i,0]+ ab*initialLocation[i,1]+ ac*initialLocation[i,2]+ optmin[3]))**2)
            func+= ((finalLocation[i,1]-(ba*initialLocation[i,0]+ bb*initialLocation[i,1]+ bc*initialLocation[i,2]+ optmin[4]))**2) 
            func+= ((finalLocation[i,2]-(ca*initialLocation[i,0]+ cb*initialLocation[i,1]+ cc*initialLocation[i,2]+ optmin[5]))**2)

#        func+=(finalLocation[i,0]-(aa*initialLocation[i,0]+ ab*initialLocation[i,1]+ ac*initialLocation[i,2]+ optmin[3]))**2 +  \
#              (finalLocation[i,1]-(ba*initialLocation[i,0]+ bb*initialLocation[i,1]+ bc*initialLocation[i,2]+ optmin[4]))**2 +  \
#              (finalLocation[i,2]-(ca*initialLocation[i,0]+ cb*initialLocation[i,1]+ cc*initialLocation[i,2]+ optmin[5]))**2


    try:
        if (rotation):
            totalObjective = math.sqrt(func)/5
        else:
            totalObjective = math.sqrt(func)/8
        print optmin
        print "totalObjective", totalObjective
        print "============================================"
        bestMatchError   = 10e16
        bestMatchError   = min([totalObjective,bestMatchError])
#        f=bestMatchError
        f=bestMatchError + 1e4*np.fabs(det - 1)
        fail = 0
        return f,g,fail        
    except Exception, e:
        print optmin,str(e)
        fail = 1
        f = 1e16
    return f,g,fail
        
cmissobject = dict()
cmissobject['numberOfWallNodes'] =  numberOfWallNodes
cmissobject['numberOfCircumfrentialNodes'] =  numberOfCircumfrentialNodes


def getSolutions(solution):
    variables = solution._variables
    if (rotation):
        consts = [0.0]*3
    else:
        consts = [0.0]*6
    for key,var in variables.iteritems():
        consts[key] = var.value
    return consts

def printObjectiveEstimate(solution,keys,filename):
    if not isinstance(solution, list):
        xv  = getSolutions(solution)
    else:
        xv = solution
    findObjective(xv,cmiss=cmissobject)
    for i,k in enumerate(keys):    
        print k,xv[i]


from pyOpt import Optimization
from pyOpt import NSGA2
from pyOpt import SLSQP
from pyOpt import ALPSO
#from pyOpt import NLPQLP    # non-linear programming sequential quadratic programming ...  
from pyOpt import MIDACO    # mized integer distributed Ant Colony Optimization
#from pyOpt import SNOPT
from pyOpt import ALHSO


if (rotation):
    x0 = [0.005, 0.005, 0.03]
else:
    x0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]


if (rotation):
    fkeys = ['rot1','rot2','rot3']
else:
    fkeys = ['rot1','rot2','rot3','translation1','translation2', 'translation3']


doGlobalSearch = False
doSwarmIntelligence = True
doMIDACO = False
doNonlinearSearch = False
doSparseNonlinear = False
doHarmonySearch = False


#print '************************************* Initial **********************************'
#printObjectiveEstimate(x0,fkeys,'BlockForInitialCondition')
#print '************************************* Initial **********************************'
if (rotation):
    bnds = [(-0.5,-0.2),(-0.3,0.3),(-0.3,0.3)]
else:
    bnds = [(-0.75,0.75),(-0.75,0.75),(-0.75,0.75),(-300.0,300.0),(-300.0,300.0), (-15.0,15.0)]

opt_prob = Optimization('GrowthOptimise',findObjective)

for i,v in enumerate(bnds):
    opt_prob.addVar(fkeys[i],'c',lower=v[0],upper=v[1],value=x0[i])
opt_prob.addObj('f')
print opt_prob

# Global Optimization
if doGlobalSearch:
    nsga2 = NSGA2()
    nsga2.setOption('maxGen', 20)
    nsga2.setOption('PopSize', 40)
    nsga2(opt_prob,cmiss=cmissobject)
    print opt_prob.solution(0)
     
    # Local Optimization Refinement
    slsqp = SLSQP()
    slsqp(opt_prob.solution(0),sens_type='FD',cmiss=cmissobject)
    print opt_prob.solution(0).solution(0)
    solution = opt_prob.solution(0).solution(0)
    
elif doSwarmIntelligence: 
    augmentedLagrnagePSO = ALPSO()
    augmentedLagrnagePSO.setOption('SwarmSize',20)
    augmentedLagrnagePSO(opt_prob,cmiss=cmissobject)      
    print opt_prob.solution(0)
     
    # Local Optimization Refinement
    slsqp = SLSQP()
    slsqp(opt_prob.solution(0),sens_type='FD',cmiss=cmissobject)
    print opt_prob.solution(0).solution(0)
    solution = opt_prob.solution(0).solution(0)
    
elif doHarmonySearch:
    augmentedLagrnageHSO = ALHSO()
    augmentedLagrnageHSO.setOption('hms',20)
    augmentedLagrnageHSO(opt_prob,cmiss=cmissobject)      
    print opt_prob.solution(0)
     
    # Local Optimization Refinement
    slsqp = SLSQP()
    slsqp(opt_prob.solution(0),sens_type='FD',cmiss=cmissobject)
    print opt_prob.solution(0).solution(0)
    solution = opt_prob.solution(0).solution(0)  

elif doSparseNonlinear: 
    sparseNonlinear = SNOPT()
    sparseNonlinear(opt_prob,cmiss=cmissobject)      
    print opt_prob.solution(0)
     
    # Local Optimization Refinement
    slsqp = SLSQP()
    slsqp(opt_prob.solution(0),sens_type='FD',cmiss=cmissobject)
    print opt_prob.solution(0).solution(0)
    solution = opt_prob.solution(0).solution(0)
    
elif doMIDACO:
    # Solve Problem (No-Parallelization)
    midaco_none = MIDACO()
    #midaco_none.setOption('IPRINT',-1)
    #midaco_none.setOption('MAXEVAL',50000)
    midaco_none(opt_prob,cmiss=cmissobject)
    print opt_prob.solution(0)
     
    # Local Optimization Refinement
    slsqp = SLSQP()
    slsqp(opt_prob.solution(0),sens_type='FD',cmiss=cmissobject)
    print opt_prob.solution(0).solution(0)
    solution = opt_prob.solution(0).solution(0)
    
   
elif doNonlinearSearch:
    # Solve Problem (No-Parallelization)
    nlpqlp_none = NLPQLP()
    nlpqlp_none.setOption('IPRINT',0)
    nlpqlp_none(opt_prob)
    print opt_prob.solution(0)
     
    # Local Optimization Refinement
    slsqp = SLSQP()
    slsqp(opt_prob.solution(0),sens_type='FD',cmiss=cmissobject)
    print opt_prob.solution(0).solution(0)
    solution = opt_prob.solution(0).solution(0)

else: 
    # Local Optimization Refinement
    slsqp = SLSQP()
    slsqp(opt_prob,sens_type='FD',cmiss=cmissobject)
    solution = opt_prob.solution(0)  
    


    
#-----------------------------------------------------------------
print '************************************* Solution **********************************'
printObjectiveEstimate(solution,fkeys,'SolutionBlock')
print '************************************* Solution **********************************'
iron.Finalise()
