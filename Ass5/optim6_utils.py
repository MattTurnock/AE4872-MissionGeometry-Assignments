#!/usr/bin/env python

"""
optim6_utils.py: Some utility functions for the OPTIM-6 assignment
Changelist:
1.1: Added convergence functionality such that the average result from the generation is used for convergence, instead
of the optimum
"""
__author__      = "Matthew Turnock"
__email__ = "matthew.turnock1@gmail.com"
__version__ = "1.1"

########################################################################################################################
import numpy as np
import pandas as pd

########################################################################################################################
#********* Misc Functions**************************

def himmelblau2(x1x2List):
    """
    Define himmelblau for use in the GA
    :param x1x2List:  a length-2 list containing x1 and x2 to be passed to Himmelblau
    :return y: The evaluated Himmelblau function
    """
    x1, x2 = x1x2List
    y = (x1**2 + x2 - 11)**2 + (x1 + x2**2 - 7)**2
    return y

#######################################################################################################################
#********** Genetic Algorithm Functions************

def getParameterData(parameterDataTemp, columns=["bitLength", "LB", "UB"]):
    """
    Function to generate standard pandas  dataframe from numpy parameter data array
    :param parameterDataTemp: A numpy-style version of the dataframe, nominally containing bit-length, LB and UB columns
    :param columns: Column labels. Nominally as shown
    :return: the Pandas dataframe representation of parameter data
    """
    parameterData = pd.DataFrame(data=parameterDataTemp, columns=columns)
    return parameterData

def get_randChrom(parameterData):
    """
    Function to create a random bit representation (chromosome), based on inputted parameterData
    :param parameterData: The pandas dataframe representation of parameter data
    :return chrom: The randomly generated chromosome
    """
    chrom = np.random.randint(2, size=sum(parameterData.loc[:, "bitLength"]))
    return chrom

def buildInitialGenome(parameterData, populationSize=100):
    """
    Function to build the initial genome (ie population)
    :param parameterData: The pandas datafram representation of parameter data
    :param populationSize: The size of population. ie number of chromosomes to generate. Default is 100
    :return genome: A randomly generated genome
    """
    chromLength = sum(parameterData.loc[:, "bitLength"])
    genome = np.zeros((populationSize, chromLength))
    for i in range(populationSize):
        genome[i, :] = get_randChrom(parameterData)
    return genome

def add_xns(parameterData, chrom):
    """
    Function to add xns to the end column of an initial parameterData dataframe
    :param parameterData: The pandas dataframe representation of parameter data
    :param chrom: The chromosome to generate the data frame from
    :return parameterDataNew: Like parameterData with an appended column for parameter values
    """
    xns = []
    for n in range(len(parameterData)):
        # Define SN1
        SN1 = sum(parameterData.loc[:n, "bitLength"])
        # Define SN2
        if n == 0:
            SN2 = 0
        else:
            SN2 = sum(parameterData.loc[:(n-1), "bitLength"])

        # Define Nn, LBn, UBn, and use them to define Delta_n
        Nn = parameterData.loc[n, "bitLength"]
        LBn = parameterData.loc[n, "LB"]
        UBn = parameterData.loc[n, "UB"]
        Delta_n = (UBn - LBn)/(2**float(Nn) - 1)

        # Define bitSum with a for loop. Note the slightly odd i behaviour from python indexing
        bitSum=0
        for i in range(SN2, SN1):
            bi = chrom[i]
            two = 2**(SN1 - (i+1))
            bitSum += (bi * two)

        # Append xn to list of xns
        xn = LBn + Delta_n*bitSum
        xns.append(xn)

    # Create new dataframe with xn appended as another column, and return it
    parameterDataNew = parameterData.copy()
    parameterDataNew.loc[:, "xn"] = xns
    return parameterDataNew

def get_ParentPairs(parents):
    """
    Function to get a list of parent pairs, given an array of parents
    :param parents: List of all parent chromosomes
    :return parentPairs: A list of randomly generated pairs of parents
    """
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

def get_children(parentPairs, splitIndex_base="random"):
    """
    Function to get the children given the list of parent pairs. Set split index to an integer if it must be enforced
    :param parentPairs: List of parent pairs to generate children from
    :param splitIndex_base: Option to use random crossover point, or impose it
    :return children: A list of the children. Same format as a genome
    """
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
    children = np.array(children)
    return children

def flipBit(bitValue):
    """
    Simple function to flip a bit
    :param bitValue: The input bit value (1 or 0)
    :return newBitValue: The flipped bit value
    """
    if bitValue == 1:
        newBitValue = 0
    elif bitValue == 0:
        newBitValue = 1
    else:
        print("Bit flip not possible, input is neither 0 or 1")

    return newBitValue


def mutateChildren(children, mutationRate=0.1/100):
    """
    Function definition to mutate the children (harsh I know but it has to be done)
    :param children: The children genome that are to be mutated (can also input any genome for mutation)
    :param mutationRate: The rate of mutation per bit
    :return mutatedChildren: The genome-like list of mutated children
    """
    # Set up the mutatedChildren array from the original one, and run a for loop through each child
    mutatedChildren = np.copy(children)
    for childIndex in range(len(children)):
        child = children[childIndex]
        # Run another for loop for each bit of the chromosome of the child
        for bitIndex in range(len(child)):
            # If the random number is less than the mutation rate, mutate the child's bit and append to the
            # mutatedChildren array
            chance = np.random.random()
            if chance < mutationRate:
                bit = child[bitIndex]
                newBit = flipBit(bit)
                mutatedChildren[childIndex][bitIndex] = newBit

    return mutatedChildren


#
def checkFitness(inputFunction, parameterDataNew, optimisation="min"):
    """
    Function definition to calculate fitness (+ve better) given an input function, parameterData including xns, and the
    optimisation mode. Note that it does not account for inputFunction that can have negative outcomes
    :param inputFunction: The input, ie objective, function that is used by the GA
    :param parameterDataNew: The pandas dataframe containing a chromosome's parameter values
    :param optimisation: The option to choose if function should be minimised or maximised
    :return fitness: The fitness of the chromosome. Bigger is always more fit
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
                       convergenceConditions=[5], printingFull=False, printRunIndex=False,
                       populationSize=100, mutationRate=0.1 / 100):
    """

    :param parameterData: Parameter data input as a pandas dataframe
    :param inputFunction: The input function given without parameters (eg himmelblau2)
    :param optimisation: The optimisation type for the given function. Eg we may want to minimise or maximise the input
    function
    :param convergenceType: Currently allowed types are "runCount", "userDefined", and "userDefinedAverage". runCount
    ends after a certain number of runs, userDefined implements a % based convergence condition, based on most optimum
    per generation; also with a minimum number of runs implemented. userDefinedAverage works similarly to the
    userDefined case, but uses the average of the previous generation instead of the optimum as convergence condition.
    :param convergenceConditions: List of convergence conditions depending on convergenceType. For "runCount" input a
    length-1 list of how many runs to complete. For "userDefined" and "userDefinedAverage" input a length-2 list with
    convergence threshold followed by minimum number of runs. TODO: implement "userDefinedAverage" and test
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
    genomeFitnessHistory = []
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
            genomeFitnessData.append(fitness)

            # Check if fitness is better than the current optimum, and update if it is
            if fitness > optimalFitness:
                optimalFitness = fitness
                optimalSolution = (chromosome, parameterDataNew, fitness)
                # Extra printing info if needed
                if printingFull:
                    print(optimalFitness)

        # Add current fitness data to the history
        genomeFitnessHistory.append(genomeFitnessData)
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
                # Define optimum of current and previous generations
                previousFitness = optimalChromosomeHistory[-2][2]
                currentFitness = optimalChromosomeHistory[-1][2]
                differenceFraction = abs((currentFitness - previousFitness) / previousFitness)
                if differenceFraction < convergenceConditions[0]:
                    convergence = True
                else:
                    convergence = False

        elif convergenceType == "userDefinedAverage":
            if (len(genomeFitnessHistory) < convergenceConditions[1]) or (len(genomeFitnessHistory) < 2):
                convergence = False

            else:
                # Define average of current and previous generations
                previousFitnessAverage = np.mean(genomeFitnessHistory[-2])
                currentFitnessAverage = np.mean(genomeFitnessHistory[-1])
                differenceFraction = abs((currentFitnessAverage - previousFitnessAverage) / previousFitnessAverage)
                if differenceFraction < convergenceConditions[0]:
                    convergence = True
                else:
                    convergence = False

                if printRunIndex:
                    differencePercent = differenceFraction*100
                    print("Difference percent: %s" %differencePercent)


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