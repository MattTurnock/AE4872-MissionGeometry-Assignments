#!/usr/bin/env python

"""test.py: Runs the Himmelblau GA for the OPTIM-6 assignment"""
__author__      = "Matthew Turnock"
__email__ = "matthew.turnock1@gmail.com"
__version__ = "1.0"

########################################################################################################################
from matplotlib import pyplot as plt
import os
import sys
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

########################################################################################################################

# A bunch of parameters for saving the files (not related to GA)
dataDirBase = "himmelblau_data_%s"
parameterDataDir = "parameterDataCSVs"
plotSavenameBase = "geneticAlgorithm_con%s_rep%s"
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
fontsz=20

# Parameters we can play with a little
convergenceConditionsBoth =  [[100], [0.1/100, 10]]
populationSize =             1000

# Often changing parameters for debugging, or changing run mode
printRunIndex =              False
printingFull =               False
show = False
save = True
titling = False
calculating = False
plotting = False
saveOutput = True

if saveOutput:
    orig_stdout = sys.stdout
    f = open('printout.txt', 'w')
    sys.stdout = f

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
            # Define directory for data to be placed
            dataDir = dataDirBase % j
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
                saveData(optimalChromosomeHistory, bitLength, convergenceTypeIndex, baseSavename=baseSavename, dataDir=dataDir,
                         parameterDataDir=parameterDataDir, parameterDataBaseSavename=parameterDataBaseSavename)

            # Load from file if calculations not performed
            else:
                optimalChromosomes, optimalFitnesses = loadData(bitLength, convergenceTypeIndex,
                                                                      baseSavename=baseSavename, dataDir=dataDir)

            # Plotting of results
            runIndices = range(len(optimalFitnesses))
            optimalValues = 1/optimalFitnesses

            # Defines which figure to plot to based on convergence type and repeat run
            if convergenceTypeIndex == 0:
                if j == 0:
                    figNumber = 1
                elif j == 1:
                    figNumber = 2
                elif j ==2:
                    figNumber = 3
            elif convergenceTypeIndex == 1:
                if j == 0:
                    figNumber = 4
                elif j == 1:
                    figNumber = 5
                elif j ==2:
                    figNumber = 6
            else: print("IT BROKE")

            # Actually does the plotting
            if plotting:
                plt.figure(figNumber)
                plt.plot(runIndices, optimalValues)
                plt.minorticks_on()
                plt.grid(which='both')
                plt.xlabel("Generation")
                plt.ylabel("Himmelblau Function Evaluation")

                legendBase = "bitLength = %s"
                plt.legend([legendBase %5, legendBase %7, legendBase %9], prop={'size': fontsz})

                ax = plt.gca()
                for item in ([ax.title, ax.xaxis.label, ax.yaxis.label, ax.yaxis.get_offset_text()] +
                             ax.get_xticklabels() + ax.get_yticklabels()):
                    item.set_fontsize(fontsz)

                saveName = plotSavenameBase %(convergenceTypeIndex, j)
                if titling: plt.title(saveName)
                if save: plt.savefig(os.path.join(plotDir, saveName  + '.pdf'),bbox_inches="tight")

            finalOptimalValue = optimalValues[-1]
            finalOptimalParameterData = add_xns(parameterData, optimalChromosomes[-1])
            finalOptimalxns = list(finalOptimalParameterData.loc[:]["xn"])

            print("\n===========================================================================")
            print("Repeat Number = %s ; Running bitLength = %s ; convergenceType = %s" % (
            j, bitLength, convergenceType))
            print("Final Himmelblau evaluation: %s" %finalOptimalValue)
            print("With parameters : \nx1 = %s \nx2 = %s" %(finalOptimalxns[0], finalOptimalxns[1]))
if show: plt.show()
# Save printed terminal output to file
if saveOutput:
    sys.stdout = orig_stdout
    f.close()