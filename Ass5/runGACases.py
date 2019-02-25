#!/usr/bin/env python

"""
runGACases .py: Runs the Himmelblau GA for the OPTIM-6 assignment
Changelist:
1.1: Implemented new average convergence option, also increased data handling efficiency TODO: explain how
"""
__author__      = "Matthew Turnock"
__email__ = "matthew.turnock1@gmail.com"
__version__ = "1.1"

########################################################################################################################
from matplotlib import pyplot as plt
from matplotlib.ticker import ScalarFormatter
import time
import matplotlib.ticker as ticker
import os
import sys
from Ass5.optim6_utils import *

def saveData(optimalChromosomeHistory, bitLength, convergenceTypeIndex, compTime, baseSavename="optimal%s_bit%s_con%s.npy",
             dataDir="", parameterDataDir="",
             parameterDataBaseSavename="optimalParameterData%s_bit%s_con%s.csv", saveLite=False):
    """
    Function that saves computed data to file - prevents need to actually run the GA many times when plotting
    :param optimalChromosomeHistory: History of optimal chromosomes provided by GA
    :param bitLength: The bit-length of the current case to be saved
    :param convergenceTypeIndex: The index of the used convergence type (only relevant in this script)
    :param compTime: The computation time for the GA which should be saved
    :param baseSavename: The base-level format for the savename of each numpy array to save (is completed later)
    :param dataDir: The directory in which to save
    :param parameterDataDir: The place to save CSV files of panda dataframes, if those are to be saved
    :param parameterDataBaseSavename: The base savename for CSV files
    :param saveLite: Boolean to determine if CSVs should be saved (implemented as there are many CSVs and it is
    time-consuming
    :return: None
    """
    optimalChromosomes, optimalParameterDatas, optimalFitnesses = optimalChromosomeHistory

    # Define savenames
    chromosomeSavename = baseSavename % ("Chromosomes", bitLength, convergenceTypeIndex)
    fitnessesSavename = baseSavename % ("Fitnesses", bitLength, convergenceTypeIndex)
    compTimeSavename = baseSavename % ("CompTime", bitLength, convergenceTypeIndex)

    # Use numpy to save the arrays
    np.save(os.path.join(dataDir, chromosomeSavename), optimalChromosomes)
    np.save(os.path.join(dataDir, fitnessesSavename), optimalFitnesses)
    np.save(os.path.join(dataDir, compTimeSavename), compTime)

    if not saveLite:
        # Save each dataframe individually
        for i in range(len(optimalParameterDatas)):
            optimalParameterData = optimalParameterDatas[i]
            parameterDataSavename = parameterDataBaseSavename % (i, bitLength, convergenceTypeIndex)
            optimalParameterData.to_csv(os.path.join(dataDir, parameterDataDir, parameterDataSavename), index=False,
                                        encoding='utf8')

# Function to load saved data from file TODO: make it load the CSV parameterData
def loadData(bitLength, convergenceTypeIndex, baseSavename="optimal%s_bit%s_con%s.npy",
             dataDir=""):
    """
    Function to load the saved data from file
    :param bitLength: The bit-length of the current case to be loaded
    :param convergenceTypeIndex: The index of the used convergence type (only relevant in this script)
    :param baseSavename: The base-level format for the savename of each numpy array to load (is completed later)
    :param dataDir: The directory in which to save
    :return optimalChromosomes, optimalFitnesses, compTime: Tuple of useful values for plotting and printing
    """

    chromosomeSavename = baseSavename % ("Chromosomes", bitLength, convergenceTypeIndex)
    fitnessesSavename = baseSavename % ("Fitnesses", bitLength, convergenceTypeIndex)
    compTimeSavename = baseSavename % ("CompTime", bitLength, convergenceTypeIndex)

    optimalChromosomes = np.load(os.path.join(dataDir, chromosomeSavename))
    optimalFitnesses = np.load(os.path.join(dataDir, fitnessesSavename))
    compTime = np.load(os.path.join(dataDir, compTimeSavename))

    return optimalChromosomes, optimalFitnesses, compTime

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
mutationRate =               0.1/100
LBBoth =                     0
UBBoth =                     5
inputFunction =              himmelblau2
optimisation =               "min"
convergenceTypes =           ["runCount", "userDefinedAverage"] #TODO: change this for new convergence conditions
fontsz=20

# Parameters we can play with a little
convergenceConditionsBoth =  [[100], [0.1/100, 5]]
populationSize =             1000

# Often changing parameters for debugging, or changing run mode
printRunIndex =              True
printingFull =               False
show =                       True
save =                       False
titling =                    False
calculating =                False
plotting =                   True
saveOutput =                 False
saveLite =                   True

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
                t0 = time.time()
                optimalChromosomeHistory = doGeneticAlgorithm(parameterData, inputFunction, optimisation=optimisation,
                                                              convergenceType=convergenceType,
                                                              convergenceConditions=convergenceConditions,
                                                              printingFull=printingFull, printRunIndex=printRunIndex,
                                                              populationSize=populationSize, mutationRate=mutationRate)
                compTime = np.array(time.time() - t0)
                # Define lists of useful data, and save to file for later loading
                optimalChromosomes, optimalParameterDatas, optimalFitnesses = optimalChromosomeHistory
                saveData(optimalChromosomeHistory, bitLength, convergenceTypeIndex, compTime, baseSavename=baseSavename,
                         dataDir=dataDir, parameterDataDir=parameterDataDir,
                         parameterDataBaseSavename=parameterDataBaseSavename, saveLite=saveLite)

            # Load from file if calculations not performed
            else:
                optimalChromosomes, optimalFitnesses, compTime = loadData(bitLength, convergenceTypeIndex,
                                                                      baseSavename=baseSavename, dataDir=dataDir)

            # Plotting of results TODO: Add plotting of number of function evaluations (indicates computation time)

            runIndices = range(len(optimalFitnesses))
            optimalValues = 1/optimalFitnesses

            generations = len(optimalFitnesses)
            computations = generations*populationSize
            # print("THIS IS IT %s" %computations)

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
                ax = plt.gca()
                ax.set_yscale('log')
                plt.minorticks_on()
                plt.grid(which='both')
                plt.xlabel("Generation")
                plt.ylabel("Himmelblau Function Evaluation")

                legendBase = "bitLength = %s"
                plt.legend([legendBase %5, legendBase %7, legendBase %9], prop={'size': fontsz})


                for item in ([ax.title, ax.xaxis.label, ax.yaxis.label, ax.yaxis.get_offset_text()] +
                             ax.get_xticklabels() + ax.get_yticklabels()):
                    item.set_fontsize(fontsz)
                ax.yaxis.set_major_formatter(ticker.FuncFormatter(
                    lambda y, pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y), 0)))).format(y)))

                saveName = plotSavenameBase %(convergenceTypeIndex, j)
                if titling: plt.title(saveName)
                if save: plt.savefig(os.path.join(plotDir, saveName  + '.pdf'),bbox_inches="tight")

                # This is where computation plots can be done - how to find out computations?

            finalOptimalValue = optimalValues[-1]
            finalOptimalParameterData = add_xns(parameterData, optimalChromosomes[-1])
            finalOptimalxns = list(finalOptimalParameterData.loc[:]["xn"])

            print("\n===========================================================================")
            print("Repeat Number = %s ; Running bitLength = %s ; convergenceType = %s" % (
            j, bitLength, convergenceType))
            print("Number of generations: %s" %len(runIndices))
            print("Computation time: %s min" %(compTime/60))
            print("Final Himmelblau evaluation: %s" %finalOptimalValue)
            print("With parameters : \nx1 = %s \nx2 = %s" %(finalOptimalxns[0], finalOptimalxns[1]))
if show: plt.show()

# Save printed terminal output to file
if saveOutput:
    sys.stdout = orig_stdout
    f.close()