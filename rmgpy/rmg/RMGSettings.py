#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2010 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
This module contains settings classes for manipulation of RMG run parameters
"""
import numpy
import logging
from rmgpy import settings

class DatabaseSettings:
    """
    class for controlling what database information is loaded
    """
    def __init__(self,
             thermoLibraries = None,
             transportLibraries = None,
             reactionLibraries = None,
             frequenciesLibraries = None,
             seedMechanisms = None,
             kineticsFamilies = 'default',
             kineticsDepositories = 'default',
             kineticsEstimator = 'rate rules',
             ):
        # This function just stores the information about the database to be loaded
        # We don't actually load the database until after we're finished reading
        # the input file
        if isinstance(thermoLibraries, str): thermoLibraries = [thermoLibraries]
        if isinstance(transportLibraries, str): transportLibraries = [transportLibraries]
        if isinstance(reactionLibraries, str): reactionLibraries = [reactionLibraries]
        if isinstance(seedMechanisms, str): seedMechanisms = [seedMechanisms]
        if isinstance(frequenciesLibraries, str): frequenciesLibraries = [frequenciesLibraries]
        self.databaseDirectory = settings['database.directory']
        self.thermoLibraries = thermoLibraries or []
        self.transportLibraries = transportLibraries
        # Modify reactionLibraries if the user didn't specify tuple input
        if reactionLibraries:
            index = 0
            while index < len(reactionLibraries):
                if isinstance(reactionLibraries[index],tuple):
                    pass
                elif isinstance(reactionLibraries[index],str):
                    reactionLibraries[index] = (reactionLibraries[index], False)
                else:
                    raise TypeError('reaction libraries must be input as tuples or strings')
                index += 1
                
        self.reactionLibraries = reactionLibraries or []
        self.seedMechanisms = seedMechanisms or []
        self.statmechLibraries = frequenciesLibraries or []
        self.kineticsEstimator = kineticsEstimator
        if kineticsDepositories == 'default':
            self.kineticsDepositories = ['training']
        elif kineticsDepositories == 'all':
            self.kineticsDepositories = None
        else:
            assert isinstance(kineticsDepositories,list), "kineticsDepositories should be either 'default', 'all', or a list of names eg. ['training','PrIMe']."
            self.kineticsDepositories = kineticsDepositories
    
        if kineticsFamilies in ('default', 'all', 'none'):
            self.kineticsFamilies = kineticsFamilies
        else:
            assert isinstance(kineticsFamilies,list), "kineticsFamilies should be either 'default', 'all', 'none', or a list of names eg. ['H_Abstraction','R_Recombination'] or ['!Intra_Disproportionation']."
            self.kineticsFamilies = kineticsFamilies
        
class ModelSettings:
    """
    class for holding the parameters affecting an RMG run
    """
    def __init__(self,toleranceMoveToCore=numpy.inf, toleranceMoveEdgeReactionToCore=numpy.inf,toleranceKeepInEdge=0.0, toleranceInterruptSimulation=1.0, 
          toleranceMoveEdgeReactionToSurface=numpy.inf, toleranceMoveSurfaceSpeciesToCore=numpy.inf, toleranceMoveSurfaceReactionToCore=numpy.inf,
          toleranceMoveEdgeReactionToSurfaceInterrupt=numpy.inf,toleranceMoveEdgeReactionToCoreInterrupt=numpy.inf, maximumEdgeSpecies=1000000, minCoreSizeForPrune=50, 
          minSpeciesExistIterationsForPrune=2, filterReactions=False, ignoreOverallFluxCriterion=False, maxNumSpecies=numpy.inf, maxNumObjsPerIter=1,terminateAtMaxObjects=False):
        
        self.fluxToleranceKeepInEdge = toleranceKeepInEdge
        self.fluxToleranceMoveToCore = toleranceMoveToCore
        self.toleranceMoveEdgeReactionToCore = toleranceMoveEdgeReactionToCore
        self.fluxToleranceInterrupt = toleranceInterruptSimulation
        self.maximumEdgeSpecies = maximumEdgeSpecies
        self.minCoreSizeForPrune = minCoreSizeForPrune
        self.minSpeciesExistIterationsForPrune = minSpeciesExistIterationsForPrune
        self.filterReactions = filterReactions
        self.ignoreOverallFluxCriterion=ignoreOverallFluxCriterion
        self.toleranceMoveEdgeReactionToSurface = toleranceMoveEdgeReactionToSurface
        self.toleranceMoveSurfaceSpeciesToCore = toleranceMoveSurfaceSpeciesToCore
        self.toleranceMoveSurfaceReactionToCore = toleranceMoveSurfaceReactionToCore
        self.terminateAtMaxObjects = terminateAtMaxObjects
        self.fluxToleranceInterrupt = toleranceInterruptSimulation
        self.toleranceMoveEdgeReactionToSurfaceInterrupt = toleranceMoveEdgeReactionToSurfaceInterrupt
        self.toleranceMoveEdgeReactionToCoreInterrupt = toleranceMoveEdgeReactionToCoreInterrupt
        self.maxNumSpecies = maxNumSpecies
        
        if maxNumObjsPerIter <= 0: #negative value or 0 value set to infinity
            logging.info('maxNumObjsPerIter was 0 or negative ... setting value to infinity')
            self.maxNumObjsPerIter = numpy.inf
        else:
            self.maxNumObjsPerIter = maxNumObjsPerIter
            
class SimulatorSettings:
    """
    class for holding the parameters affecting the behavior of the solver
    """
    def __init__(self,atol=1e-16, rtol=1e-8, sens_atol=1e-6, sens_rtol=1e-4):
        self.atol = atol
        self.rtol = rtol
        self.sens_atol = sens_atol
        self.sens_rtol = sens_rtol