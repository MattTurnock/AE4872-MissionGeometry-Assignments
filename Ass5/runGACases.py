import numpy as np
from matplotlib import pyplot as plt
import os
from Ass5.optim6_utils import *

# Function to save data to file (because the calculations take FOREVER)
def saveData(optimalChromosomeHistory, bitLength, convergenceTypeIndex, baseSavename="optimal%s_bit%s_con%s.npy",
             dataDir="", parameterDataDir="",
             parameterDataBaseSavename="optimalParameterData%s_bit%s_con%s.csv"):
    optimalChromosomes, optimalParameterDatas, optimalFitnesses = optimalChromosomeHistory

    # Define savenames
    chromosomeSavename = baseSavename % ("Chromosomes", bitLength, convergenceTypeIndex)
    fitnessesSavename = baseSavename % ("Fitnesses", bitLength, convergenceTypeIndex)

    # Use numpy to save the arrays
    np.save(os.path.join(dataDir, chromosomeSavename), optimalChromosomes)
    np.save(os.path.join(dataDir, fitnessesSavename), optimalFitnesses)

    # Save each dataframe individually
    for i in range(len(optimalParameterDatas)):
        optimalParameterData = optimalParameterDatas[i]
        parameterDataSavename = parameterDataBaseSavename % (i, bitLength, convergenceTypeIndex)
        optimalParameterData.to_csv(os.path.join(dataDir, parameterDataDir, parameterDataSavename), index=False,
                                    encoding='utf8')

# Function to load saved data from file TODO: make it load the CSV parameterData
def loadData(bitLength, convergenceTypeIndex, baseSavename="optimal%s_bit%s_con%s.npy",
             dataDir=""):

    chromosomeSavename = baseSavename % ("Chromosomes", bitLength, convergenceTypeIndex)
    fitnessesSavename = baseSavename % ("Fitnesses", bitLength, convergenceTypeIndex)

    optimalChromosomes = np.load(os.path.join(dataDir, chromosomeSavename))
    optimalFitnesses = np.load(os.path.join(dataDir, fitnessesSavename))

    return optimalChromosomes, optimalFitnesses



# A bunch of parameters for saving the files (not related to GA)
dataDirBase = "himmelblau_data_%s"
parameterDataDir = "parameterDataCSVs"
plotDir = "plots"
CSVDir = "parameterDataCSVs"
baseSavename = "optimal%s_bit%s_con%s.npy"
parameterDataBaseSavename = "optimalParameterData%s_bit%s_con%s.csv"

# Parameters that are not to be changed (are hard given in assignment)
bitLengths =                 [5, 7, 9]
repeats =                    3
mutationRate =               0.1/1000
LBBoth =                     0
UBBoth =                     5
inputFunction =              himmelblau2
optimisation =               "min"
convergenceTypes =           ["runCount", "userDefined"]

# Parameters we can play with a little
convergenceConditionsBoth =  [[100], [0.1/100, 10]]
populationSize =             1000

# Often changing parameters for debugging, or changing run mode
printRunIndex =              True
printingFull =               False
calculating = True

# Loop to repeat the GA repeats times
for j in range(repeats):

    # Loop to do each bit-length
    for bitLength in bitLengths:

        # Some print statements to let the user know what type of run is being completed
        if calculating:
            toprint = "calculation"
        else:
            toprint = "load from file"
        print("Method is: %s \n" %toprint)

        # Initialise the parameterData with a numpy array, then generate the dataFrame
        parameterDataTemp = np.array([[bitLength, LBBoth, UBBoth],
                                       [bitLength, LBBoth, UBBoth]])

        parameterData = getParameterData(parameterDataTemp)

        # For loop to run through both types of convergence (ie forced and 'true' convergence)
        for convergenceTypeIndex in range(len(convergenceTypes)):
            # define the convergence type and print a statement saying what type of run is being done
            convergenceType = convergenceTypes[convergenceTypeIndex]
            print("Repeat Number = %s ; Running bitLength = %s ; convergenceType = %s\n" % (j, bitLength, convergenceType))

            # Only actually perform the GA if calculating is true, otherwise load from file
            if calculating:

                # Define convergence conditions and run the GA
                convergenceConditions = convergenceConditionsBoth[convergenceTypeIndex]
                optimalChromosomeHistory = doGeneticAlgorithm(parameterData, inputFunction, optimisation=optimisation,
                                                              convergenceType=convergenceType,
                                                              convergenceConditions=convergenceConditions,
                                                              printingFull=printingFull, printRunIndex=printRunIndex,
                                                              populationSize=populationSize, mutationRate=mutationRate)

                # Define lists of useful data, and save to file for later loading
                optimalChromosomes, optimalParameterDatas, optimalFitnesses = optimalChromosomeHistory
                dataDir = dataDirBase %j
                saveData(optimalChromosomeHistory, bitLength, convergenceTypeIndex, baseSavename=baseSavename, dataDir=dataDir,
                         parameterDataDir=parameterDataDir, parameterDataBaseSavename=parameterDataBaseSavename)

            # Load from file if calculations not performed
            else:
                optimalChromosomes, optimalFitnesses = loadData(bitLength, convergenceTypeIndex,
                                                                      baseSavename=baseSavename, dataDir=dataDir)

            # TODO:Do some plotting
            runIndices = range(len(optimalFitnesses))
            print(runIndices)
            print(optimalFitnesses)

            # plt.figure()
            # plt.plot(runIndices, optimalFitnesses)
            # plt.show()

            # print("PROPER : %s\n" %optimalChromosomes)
            # print("NEW    : %s\n" %optimalChromosomesNew)

# print("REGULAR : \n%s" %optimalParameterDatas[-1])
# print("LOADED : \n%s" %pd.read_csv("output_filename.csv"))
