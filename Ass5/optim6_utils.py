#!/usr/bin/env python

"""test.py: Some utility functions for the OPTIM-6 assignment"""
__author__      = "Matthew Turnock"
__email__ = "matthew.turnock1@gmail.com"
__version__ = "1.0"

########################################################################################################################
import numpy as np
import pandas as pd
import itertools, random
from Ass4.optim1_utils import himmelblau

########################################################################################################################
#********* Misc Functions**************************
# Define himmelblau for use in the GA
def himmelblau2(x1x2List):
    x1, x2 = x1x2List
    return (x1**2 + x2 - 11)**2 + (x1 + x2**2 - 7)**2


#######################################################################################################################
#********** Genetic Algorithm Functions************

# Function to generate standard  dataframe from numpy parameter data array
def getParameterData(parameterDataTemp, columns=["bitLength", "LB", "UB"]):
    parameterData = pd.DataFrame(data=parameterDataTemp, columns=columns)
    return parameterData

# Function to create a random bit representation (chromosome), based on inputted parameterData
def get_randChrom(parameterData):
    chrom = np.random.randint(2, size=sum(parameterData.loc[:, "bitLength"]))
    return chrom

def buildInitialGenome(parameterData, populationSize=100):
    chromLength = sum(parameterData.loc[:, "bitLength"])
    genome = np.zeros((populationSize, chromLength))
    for i in range(populationSize):
        genome[i, :] = get_randChrom(parameterData)
    return genome

# Function to add xns to the end column of an initial parameterData dataframe
def add_xns(parameterData, chrom):

    xns = []
    for n in range(len(parameterData)):
        # Define SN1
        SN1 = sum(parameterData.loc[:n, "bitLength"])

        # Define SN2
        if n == 0:
            SN2 = 0
        else:
            SN2 = sum(parameterData.loc[:(n-1), "bitLength"])
            # print("SN2 is %s" % SN2)

        # Define Nn, LBn, UBn, and use them to define Delta_n
        Nn = parameterData.loc[n, "bitLength"]
        LBn = parameterData.loc[n, "LB"]
        UBn = parameterData.loc[n, "UB"]
        Delta_n = (UBn - LBn)/(2**float(Nn) - 1)
        # print("IT IS %s" %Delta_n)

        # Define bitSum with a for loop. Note the slightly odd i behaviour from python indexing
        bitSum=0
        for i in range(SN2, SN1):
            bi = chrom[i]
            # print("Nn is: %s. i is: %s. power is %s" %(Nn, i, SN1 - (i+1)))
            two = 2**(SN1 - (i+1))
            bitSum += (bi * two)

        # Append xn to list of xns
        xn = LBn + Delta_n*bitSum
        xns.append(xn)

    # Create new dataframe with xn appended as another column, and return it
    parameterDataNew = parameterData.copy()
    parameterDataNew.loc[:, "xn"] = xns
    return parameterDataNew

# Function to get a list of parent pairs, given an array of parents
def get_ParentPairs(parents):

    # If there are an odd number of parents, one cannnot make pairs from them, so dont allow it
    if len(parents)%2 != 0:
        print("Length of Parents is not even, illegal operation")

    # Make list of parent indices for later, and an empty list for the parent pairs
    parentIndices = list(range(len(parents)))
    parentPairs = []

    # While loop until parentIndices is empty
    running = True
    while running:
        # Choose first parent, then delete index. Do same for second parent. Add both parents to pairs list
        parent1 = np.random.choice(parentIndices)
        parentIndices.remove(parent1)
        parent2 = np.random.choice(parentIndices)
        parentIndices.remove(parent2)
        parentPairs.append((parents[parent1], parents[parent2]))

        # When parentIndices is empty, loop is over
        if len(parentIndices) == 0: running = False

    return parentPairs

# Function to get the children given the list of parent pairs. Set split index to an integer if it must be enforced
def get_children(parentPairs, splitIndex_base="random"):

    # Open empty array to dump the new children
    children = []
    # For loop that iterates through parentPairs and mixes them according to splitIndex
    for i in range(len(parentPairs)):
        # Create pair and define parents
        pair = parentPairs[i]
        parent1, parent2 = pair

        # Check for randomly generated splitIndex or enforced
        if splitIndex_base == "random":
            # List all possible split indices
            splitIndices = list(range(len(parent1)))
            # Delete first and end entries (would result in children identical to parents)
            del splitIndices[0]; del splitIndices[-1]
            # Randomly pick splitIndex from the list
            splitIndex = np.random.choice(splitIndices)
        else: splitIndex = splitIndex_base

        # Define first and second halves of each parent's chromosomes
        parent1_first, parent1_second = np.split(parent1, [splitIndex])
        parent2_first, parent2_second = np.split(parent2, [splitIndex])

        # Make the children and append to list
        child1 = np.concatenate((parent1_first, parent2_second))
        child2 = np.concatenate((parent2_first, parent1_second))
        children.append(child1)
        children.append(child2)
    # Return children as an array
    return np.array(children)

# Simple function to flip a bit
def flipBit(bitValue):
    if bitValue == 1:
        newBitValue = 0
    elif bitValue == 0:
        newBitValue = 1
    else:
        print("Bit flip not possible, input is neither 0 or 1")

    return newBitValue


# Function definition to mutate the children (harsh I know but it has to be done)
def mutateChildren(children, mutationRate=0.1/100):
    # Set up the mutatedChildren array from the original one, and run a for loop through each child
    mutatedChildren = np.copy(children)
    for childIndex in range(len(children)):
        child = children[childIndex]
        # Run another for loop for each bit of the chromosome of the child
        for bitIndex in range(len(child)):
            # If the random number is less than the mutation rate, mutate the child's bit and append to the mutatedChildren array
            chance = np.random.random()
            if chance < mutationRate:
                bit = child[bitIndex]
                newBit = flipBit(bit)
                mutatedChildren[childIndex][bitIndex] = newBit

    return mutatedChildren


# Function definition to calculate fitness (+ve better) given an input function, parameterData including xns, and the optimisation mode
def checkFitness(inputFunction, parameterDataNew, optimisation="min"):
    """
    IMPORTANT NOTE:
    Currently this function does not account for inputFunction that can have -ve outcomes
    """

    # Put xns in a simple list to give to the input function
    xns = parameterDataNew.loc[:, "xn"]
    # Pass xns to input function. Note that input to the function MUST be as a list in its own definition also
    output = inputFunction(xns)

    # For minimisation, take inverse solution, otherwise take the regular output
    if optimisation == "min":
        fitness = 1 / output
    elif optimisation == "max":
        fitness = output

    return fitness


def doGeneticAlgorithm(parameterData, inputFunction, optimisation="min", convergenceType="runCount",
                       convergenceConditions=[5], printingFull=False, printRunIndex=False, populationSize=100, mutationRate=0.1 / 100):
    """

    :param parameterData: Parameter data input as a pandas dataframe
    :param inputFunction: The input function given without parameters (eg himmelblau2)
    :param optimisation: The optimisation type for the given function. Eg we may want to minimise or maximise the input
    function
    :param convergenceType: Currently allowed types are "runCount" and "userDefined". runCount ends after a certain
    number of runs, userDefined implements a % based convergence condition, also with a maximum number of allowed
    identical iterations (SECOND PART NOT CURRENTLY IMPLEMENTED)
    :param convergenceConditions: List of convergence conditions depending on convergenceType
    :param printingFull: Adds more printing info, usually leave false
    :param printingRunIndex: Simply prints the runIndex to give an idea of how far through we are
    :param populationSize: Size of the population
    :param mutationRate: Rate of mutation per bit of a chromosome
    :return optimalChromosomeHistory: A list of info for later use containing tuples of
    (chromosome, parameterDataNew, fitness)
    """

    # Build initial random genome
    initialGenome = buildInitialGenome(parameterData, populationSize=populationSize)

    # Intialise a bunch of parameters for use in later loops
    optimalChromosomeHistory = []
    genome = initialGenome
    optimalFitness = 0
    runIndex = 0
    convergence = False

    # While loop to run until convergence is achieved
    while not convergence:
        # Add to the runIndex for convergence cases
        runIndex += 1
        if printRunIndex:
            print("Current runIndex: %s" %runIndex)
        # Initialise empty list to dump fitness data into (may be useful later)
        genomeFitnessData = []

        # For loop to find the fitness of each chromosome in the genome, and update the most optimal solution
        for chromosome in genome:
            # Add xns to final column of parameterData
            parameterDataNew = add_xns(parameterData, chromosome)
            # Calculate fitness with predefined function
            fitness = checkFitness(inputFunction, parameterDataNew, optimisation=optimisation)
            # Add fitness data to the fitnessData list
            genomeFitnessData.append((chromosome, parameterDataNew, fitness))

            # Check if fitness is better than the current optimum, and update if it is
            if fitness > optimalFitness:
                optimalFitness = fitness
                optimalSolution = (chromosome, parameterDataNew, fitness)
                # Extra printing info if needed
                if printingFull:
                    print(optimalFitness)

        # Add optimal solution to the history list
        optimalChromosomeHistory.append(optimalSolution)
        # Use predefined function to make children and mutate them
        parentPairs = get_ParentPairs(genome)
        children = get_children(parentPairs)
        mutatedChildren = mutateChildren(children, mutationRate=mutationRate)
        # Update genome with mutated children
        genome = mutatedChildren

        if printingFull:
            print(optimalChromosomeHistory)

        # Check for runCount based convergence
        if convergenceType == "runCount":
            if runIndex > convergenceConditions[0]:
                convergence = True
            else:
                convergence = False

        # A user defined convergence condition has 2 components: convergence % and minimum number of allowed iterations
        # (which must be at least 2)
        elif convergenceType == "userDefined":
            # If the history has length less than 2, then convergence cannot have happened
            if (len(optimalChromosomeHistory) < convergenceConditions[1]) or (len(optimalChromosomeHistory) < 2):
                convergence = False
            else:
                previousFitness = optimalChromosomeHistory[-2][2]
                currentFitness = optimalChromosomeHistory[-1][2]
                differenceFraction = abs((currentFitness - previousFitness) / previousFitness)
                if differenceFraction < convergenceConditions[0]:
                    convergence = True
                else:
                    convergence = False

    optimalChromosomes = []
    optimalParameterDatas = []
    optimalFitnesses = []

    for case in optimalChromosomeHistory:
        optimalChromosomes.append(case[0])
        optimalParameterDatas.append(case[1])
        optimalFitnesses.append(case[2])

    optimalChromosomes = np.array(optimalChromosomes)
    optimalFitnesses = np.array(optimalFitnesses)

    return (optimalChromosomes, optimalParameterDatas, optimalFitnesses)


########################################################################################################################
# Testing code below
#
# np.random.seed(12)
# parameterData_temp = np.array([[3, 0, 5],
#                                [32, 0, 5],
#                                [4, 0, 5]])
#
# parameterData = pd.DataFrame(data=parameterData_temp, columns=["bitLength", "LB", "UB"])
# # print(parameterData.loc[:0, "bitLength"])
# # parameterData.loc[:,"cock"] = [1,2,3]
# print(parameterData)
# chromLength = sum(parameterData.loc[:, "bitLength"])
# # print(chromLength)
#
# chromCount = 100
# chromArray = np.zeros((chromCount, chromLength))
# for i in range(chromCount):
#
#     chromSample = get_randChrom(parameterData)
#     chromArray[i, :] = chromSample
#     # bitRepSamples = np.vstack((bitRepSamples, bitRepSample))
#     # print(add_xns(parameterData, bitRepSample))
#
# parentPairs = get_ParentPairs(chromArray)
# children = get_children(parentPairs)
#
# mutationRate = 0.1/100
# # print(children)
# mutatedChildren = mutateChildren(children, mutationRate=mutationRate)
# # print(mutatedChildren)
# #
# # print("COCK")
# # print(buildInitialGenome(parameterData, populationSize=10))
#
# testChrom = mutatedChildren[0]
# print(testChrom)
#
# def testFunction(x1x2x3List):
#     x1, x2, x3 = x1x2x3List
#     y = 2*x1**2 - 3*x1*x2 + 6*x3**3 -7
#     return abs(y)
#
# parameterDataNew = add_xns(parameterData, testChrom)
# print(checkFitness(testFunction, parameterDataNew))
# print(parameterDataNew)

# START HERE TRYING TO WRITE THE FULL FUNCTION FOR GENETIC ALGORITHMOS
# --> First initialise the population, then go ahead with fitness claculations,
#     will need to add functions for full population fitness checks and for convergence checks etc

# from matplotlib import pyplot as plt
#
# # np.random.seed(12)
# parameterData_temp = np.array([[32, 0, 5],
#                                [32, 0, 5]])
#
# parameterData = pd.DataFrame(data=parameterData_temp, columns=["bitLength", "LB", "UB"])
#
# inputFunction = himmelblau2
# optimisation = "min"
# # convergenceType = "userDefined"
# # convergenceConditions = [0.1/100, 3]  # Fraction followed by maximum identical iterations
# convergenceType = "runCount"
# convergenceConditions = [100]  # This one is used for if convergence type is by run Number
#
# printingFull = False
# printRunIndex = True
# populationSize = 100
# mutationRate = 0.1/100
#
# geneticHistory = doGeneticAlgorithm(parameterData, himmelblau2, optimisation=optimisation, convergenceType=convergenceType,
#                        convergenceConditions=convergenceConditions, printingFull=printingFull, printRunIndex=printRunIndex,
#                          populationSize=populationSize, mutationRate=mutationRate)
#
# fitnesses = []
# iterations = range(len(geneticHistory))
# for i in iterations:
#     fitnesses.append(geneticHistory[i][2])
#
# plt.figure()
# plt.plot(iterations, fitnesses)
# plt.show()