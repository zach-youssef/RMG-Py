###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""
This module contains graph isomorphism functions implementing the
VF3 algorithm of Vento and Foggia (University of Salerno, Italy)
"""

from __future__ import division
cimport cython
from rmgpy.molecule.graph import Graph
from rmgpy.exceptions import VF2Error
from rmgpy.molecule.vf2 import VF2

################################################################################

cdef class VF3(VF2):
    """
    An implementation of the third version of Vento-Foggia (VF3) algorithm
    for graph and subgraph isomorphism

    The core structure of the algorithm is the same as VF2, but with the added improvements
    introduced in VF2 Plus and VF3.

    The improvements include pre-sorting the first graph into a search order using the
    probabilities of vertex characteristics, then using this ordering to pre-determine
    the core and feasibility sets of the first graph at each level of the search state.

    It also uses an improved method to select candidate nodes from the second graph.
    """

    cdef preprocess(self, Graph graph1, Graph graph2, dict init_map, bint subgraph):
        """
        Sorts the vertices of the smaller graph into the order they will be searched in,
        then uses this ordering to pre-determine the core and feasibility sets for each level
        of the state-search.
        """
        cdef Vertex vertex1, vertex2, neighbor
        cdef dict labelCounts, degreeCounts
        cdef int graph1_size

        cdef int degree, max_degree

        cdef float p_label, p_degree, prob

        cdef dict p_feasible

        cdef list exploration_sequence

        cdef int level

        cdef list terminal_levels  #at each call depth, the list of nodes that enter the terminal set
        cdef list core_ordering  #the order the nodes from graph2 will enter the mapping

        cdef dict parent_map  # Map from each node to its parent

        graph1_size = len(graph1.vertices)

        # Count the number of times each label and degree appears in the larger graph
        for vertex1 in graph1.vertices:
            labelCounts[getattr(vertex1, 'label', "")] += 1
            degreeCounts[len(vertex1.edges)] += 1

        # Find the probability of a node with a compatible label and degree
        # appearing in graph1 for each node in graph2
        if subgraph:
            # Find the maximum degree so we know what degrees we need to consider
            # when calculating probabilities on the vertices of the subgraph
            max_degree = 0

            for degree in degreeCounts.keys():
                if degree > max_degree:
                    max_degree = degree

        for vertex2 in graph2.vertices:
            p_label = labelCounts.get(getattr(vertex2, 'label', ""), 0) / graph1_size

            if subgraph:
                p_degree = 0
                for degree in range(len(vertex2.edges), max_degree + 1):
                    p_degree += degreeCounts.get(degree, 0) / graph1_size
            else:
                p_degree = degreeCounts.get(vertex2.label, 0) / graph1_size

            p_feasible[vertex2] = p_degree * p_label

        #Create a sorted list of graph2's vertices from low -> high probability
        #Nodes within the initial mapping are moved to the front of the sequence
        cdef sort_label(Vertex vertex2):
            if init_map is not None and vertex2 in init_map.values():
                return 0
            else:
                return p_feasible[vertex2]

        exploration_sequence = sorted(graph2.vertices, key=sort_label)

        #Using the exploration sequence, pre-process the terminal sets of graph2 for each level of the matching
        level = 0

        for vertex2 in exploration_sequence:
            terminal_levels.append([])
            for neighbor in vertex2.edges:
                terminal_levels[level].append(neighbor)
            level += 1

        #Store the results of the pre-processing
        self.parents = parent_map
        self.exploration_sequence = exploration_sequence
        self.terminal_levels = terminal_levels

        return

    cdef isomorphism(self, Graph graph1, Graph graph2, dict initialMapping, bint subgraph, bint findAll):
        """
        Evaluate the isomorphism relationship between graphs `graph1` and
        `graph2` with optional initial mapping `initialMapping`. If `subgraph`
        is ``True``, `graph2` is treated as a possible subgraph of `graph1`.
        If `findAll` is ``True``, all isomorphisms are found; otherwise only
        the first is found.
        """

        #Peform pre-processing
        self.preprocess(self, graph1, graph2, initialMapping, subgraph)

        #The rest of the setup before calling match is the same
        #super.isomorphism will then call VF3's new match function
        super(VF3, self).isomorphism(self, graph1, graph2, initialMapping, subgraph, findAll)

    cdef bint match(self, int callDepth) except -2:
        """
        Recursively search for pairs of vertices to match, until all vertices
        are matched or the viable set of matches is exhausted. The `callDepth`
        parameter helps ensure we never enter an infinite loop.
        """

        cdef dict mapping
        cdef list candidates
        cdef bint has_terminal
        cdef Vertex vertex1, vertex2

        # The call depth should never be negative!
        if callDepth < 0:
            raise VF2Error('Negative call depth encountered in VF3_match().')

        # Done if we have mapped to all vertices in graph
        if callDepth == 0:
            if self.findAll:
                mapping = {}
                for vertex2 in self.graph2.vertices:
                    if vertex2.ignore:
                        continue
                    assert vertex2.mapping is not None
                    assert vertex2.mapping.mapping is vertex2
                    mapping[vertex2.mapping] = vertex2
                self.mappingList.append(mapping)
            self.isMatch = True
            return True

        #Update terminal status of graph2 vertices
        for vertex2 in self.terminal_levels[len(self.graph2.vertices) - callDepth]:
            vertex2.terminal = True

        #Next candidate node from graph2 is predetermined
        vertex2 = self.exploration_sequence[len(self.graph2.vertices) - callDepth]

        #Graph1 candidates for the pairing are the neighbors of the node mapped to the 'parent'
        # of vertex2,if possible
        candidates = self.parents[vertex2].mapping.edges if self.parents.get(vertex2, False) \
            else self.graph1.vertices

        #Store whether vertex2 is terminal so we don't have to re-calculate when we backtrack
        has_terminal = vertex2.terminal

        #Search candidate pairs in graph1
        for vertex1 in candidates:
            #Ignore if marked as ignore, or if already in mapping
            if vertex1.ignore or vertex1.mapping is not None:
                continue
            #If terminal nodes are available we only search those
            if has_terminal and not vertex1.terminal:
                continue
            #Matching nodes will have the same label
            if getattr(vertex2, 'label', '') is not getattr(vertex1, 'label', ''):
                continue
            #Propose pairing
            if self.feasible(vertex1, vertex2):
                # Add proposed match to mapping
                self.addToMapping(vertex1, vertex2)
                # Recurse
                isMatch = self.match(callDepth - 1)
                if isMatch and not self.findAll:
                    return True
                # Undo proposed match
                self.removeFromMapping(vertex1, vertex2)

        # If none of the matches lead to a complete isomorphism, return false
        return False

    cdef addToMapping(self, Vertex vertex1, Vertex vertex2):
        """
        Adds valid mapping from vertex1 to vertex2, but only updates adds terminal status
        for vertices in graph1 
        
        (graph2 terminal status at each level is determined in pre-processing and applied in match)
        """

        cdef Vertex v

        #Map the vertices to one another
        # Map the vertices to one another
        vertex1.mapping = vertex2
        vertex2.mapping = vertex1

        # Remove these vertices from the set of terminals
        vertex1.terminal = False
        vertex2.terminal = False

        #Add unmapped neighbors of vertex1 to the terminal set
        for v in vertex1.edges:
            v.terminal = v.mapping is None

    cdef removeFromMapping(self, Vertex vertex1, Vertex vertex2):
        """
        Removes a mapping between vertex2, and updates terminal status of graph1 accordingly.
        
        (graph2 terminal status is still handled by the preprocessed data)
        """

        cdef Vertex v, v2

        #Unmap the vertices
        vertex1.mapping = None
        vertex2.mapping = None

        # Recompute the terminal status of any neighboring atoms
        for v in vertex1.edges:
            if v.mapping is not None: continue
            for v2 in v.edges:
                if v2.mapping is not None:
                    v.terminal = True
                    break
            else:
                v.terminal = False

        # Restore terminal status
        for v in vertex1.edges:
            if v.mapping is not None:
                vertex1.terminal = True
                break
            else:
                vertex1.terminal = False

        # Recompute the terminal status of any neighboring atoms
        for v in vertex1.edges:
            if v.mapping is not None: continue
            for v2 in v.edges:
                if v2.mapping is not None:
                    v.terminal = True
                    break
            else:
                v.terminal = False