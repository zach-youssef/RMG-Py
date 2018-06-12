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

cdef class VF3:
    """
    An implementation of the third version of Vento-Foggia (VF3) algorithm
    for graph and subgraph isomorphism

    The core structure of the algorithm is the same as VF2, but with the added improvements
    introduced in VF2 Plus and VF3.

    The improvements include pre-sorting the first graph into a search order using the
    probabilities of vertex characteristics, then using this ordering to pre-determine
    the core and feasibility sets of the first graph at each level of the search state.

    It also uses an improved method to select candidate nodes from the second graph.

    Code within the first #### enclosed block is identical to VF2 code,
    extending VF2 to reuse the code refused to compile
    """
    ################################################################################

    def __init__(self, graphA = None, graphB = None):
        self.graph1 = graphA
        self.graph2 = graphB

    @property
    def graphA(self):
        return self.graph1

    @graphA.setter
    def graphA(self, value):
        self.graph1 = value
        self.graph1.sortVertices()

    @property
    def graphB(self):
        return self.graph2

    @graphB.setter
    def graphB(self, value):
        self.graph2 = value
        self.graph2.sortVertices()

    cpdef bint isIsomorphic(self, Graph graph1, Graph graph2, dict initialMapping) except -2:
        """
        Return ``True`` if graph `graph1` is isomorphic to graph `graph2` with
        the optional initial mapping `initialMapping`, or ``False`` otherwise.
        """
        self.isomorphism(graph1, graph2, initialMapping, False, False)
        return self.isMatch

    cpdef list findIsomorphism(self, Graph graph1, Graph graph2, dict initialMapping):
        """
        Return a list of dicts of all valid isomorphism mappings from graph
        `graph1` to graph `graph2` with the optional initial mapping 
        `initialMapping`. If no valid isomorphisms are found, an empty list is
        returned.
        """
        self.isomorphism(graph1, graph2, initialMapping, False, True)
        return self.mappingList

    cpdef bint isSubgraphIsomorphic(self, Graph graph1, Graph graph2, dict initialMapping) except -2:
        """
        Return ``True`` if graph `graph1` is subgraph isomorphic to subgraph
        `graph2` with the optional initial mapping `initialMapping`, or
        ``False`` otherwise.
        """
        self.isomorphism(graph1, graph2, initialMapping, True, False)
        return self.isMatch

    cpdef list findSubgraphIsomorphisms(self, Graph graph1, Graph graph2, dict initialMapping):
        """
        Return a list of dicts of all valid subgraph isomorphism mappings from
        graph `graph1` to subgraph `graph2` with the optional initial mapping 
        `initialMapping`. If no valid subgraph isomorphisms are found, an empty
        list is returned.
        """
        self.isomorphism(graph1, graph2, initialMapping, True, True)
        return self.mappingList

    cdef isomorphism(self, Graph graph1, Graph graph2, dict initialMapping, bint subgraph, bint findAll):
        """
        Evaluate the isomorphism relationship between graphs `graph1` and
        `graph2` with optional initial mapping `initialMapping`. If `subgraph`
        is ``True``, `graph2` is treated as a possible subgraph of `graph1`.
        If `findAll` is ``True``, all isomorphisms are found; otherwise only
        the first is found.
        """
        cdef int callDepth, index1, index2

        if self.graph1 is not graph1:
            self.graph1 = graph1
            graph1.sortVertices()

        if self.graph2 is not graph2:
            self.graph2 = graph2
            graph2.sortVertices()

        #Peform pre-processing
        self.preprocess(graph1, graph2, initialMapping, subgraph)

        self.initialMapping = initialMapping
        self.subgraph = subgraph
        self.findAll = findAll

        # Clear previous result
        self.isMatch = False
        self.mappingList = []

        # Some quick isomorphism checks based on graph sizes
        if not self.subgraph and len(graph2.vertices) != len(graph1.vertices):
            # The two graphs don't have the same number of vertices, so they
            # cannot be isomorphic
            return
        elif not self.subgraph and len(graph2.vertices) == len(graph1.vertices) == 0:
            # The two graphs don't have any vertices; this means they are
            # trivially isomorphic
            self.isMatch = True
            return
        elif self.subgraph and len(graph2.vertices) > len(graph1.vertices):
            # The second graph has more vertices than the first, so it cannot be
            # a subgraph of the first
            return

        # Initialize callDepth with the size of the smallest graph
        # Each recursive call to match() will decrease it by one;
        # when the whole graph has been explored, it should reach 0
        # It should never go below zero!
        callDepth = len(graph2.vertices)

        # Initialize mapping by clearing any previous mapping information
        for vertex1 in graph1.vertices:
            vertex1.mapping = None
            vertex1.terminal = False
        for vertex2 in graph2.vertices:
            vertex2.mapping = None
            vertex2.terminal = False
        # Set the initial mapping if provided
        if self.initialMapping is not None:
            for vertex1, vertex2 in self.initialMapping.items():
                self.addToMapping(vertex1, vertex2)
            callDepth -= len(self.initialMapping)

        self.match(callDepth)

    cpdef bint feasible(self, Vertex vertex1, Vertex vertex2) except -2:
        """
        Return ``True`` if vertex `vertex1` from the first graph is a feasible
        match for vertex `vertex2` from the second graph, or ``False`` if not.
        The semantic and structural relationship of the vertices is evaluated,
        including several structural "look-aheads" that cheaply eliminate many
        otherwise feasible pairs.
        """
        cdef Vertex vert1, vert2
        cdef Edge edge1, edge2
        cdef int term1Count, term2Count, neither1Count, neither2Count

        if not self.subgraph:
            # To be feasible the connectivity values must be an exact match
            if vertex1.connectivity1 != vertex2.connectivity1: return False
            if vertex1.connectivity2 != vertex2.connectivity2: return False
            if vertex1.connectivity3 != vertex2.connectivity3: return False

        # Semantic check #1: vertex1 and vertex2 must be equivalent
        if self.subgraph:
            if not vertex1.isSpecificCaseOf(vertex2): return False
        else:
            if not vertex1.equivalent(vertex2): return False

        # Semantic check #2: adjacent vertices to vertex1 and vertex2 that are
        # already mapped should be connected by equivalent edges
        for vert2 in vertex2.edges:
            if vert2.mapping is not None:
                vert1 = vert2.mapping
                if vert1 not in vertex1.edges:
                    # The vertices are joined in graph2, but not in graph1
                    return False
                edge1 = vertex1.edges[vert1]
                edge2 = vertex2.edges[vert2]
                if self.subgraph:
                    if not edge1.isSpecificCaseOf(edge2): return False
                else:
                    if not edge1.equivalent(edge2): return False

        # There could still be edges in graph1 that aren't in graph2; this is okay
        # for subgraph matching, but not for exact matching
        if not self.subgraph:
            for vert1 in vertex1.edges:
                if vert1.mapping is not None:
                    if vert1.mapping not in vertex2.edges:
                        # The vertices are joined in graph1, but not in graph2
                        return False

        # Count number of terminals adjacent to vertex1 and vertex2
        term1Count = 0;
        term2Count = 0;
        neither1Count = 0;
        neither2Count = 0
        for vert1 in vertex1.edges:
            if vert1.terminal:
                term1Count += 1
            elif vert1.mapping is not None:
                neither1Count += 1
        for vert2 in vertex2.edges:
            if vert2.terminal:
                term2Count += 1
            elif vert2.mapping is not None:
                neither2Count += 1

        # Level 2 look-ahead: the number of adjacent vertices of vertex1 and
        # vertex2 that are non-terminals must be equal
        if self.subgraph:
            if neither1Count < neither2Count: return False
        else:
            if neither1Count != neither2Count: return False

        # Level 1 look-ahead: the number of adjacent vertices of vertex1 and
        # vertex2 that are terminals must be equal
        if self.subgraph:
            if term1Count < term2Count: return False
        else:
            if term1Count != term2Count: return False

        # Level 0 look-ahead: all adjacent vertices of vertex2 already in the
        # mapping must map to adjacent vertices of vertex1
        for vert2 in vertex2.edges:
            if vert2.mapping is not None:
                if vert2.mapping not in vertex1.edges: return False
        # Also, all adjacent vertices of vertex1 already in the mapping must map to
        # adjacent vertices of vertex2, unless we are subgraph matching
        if not self.subgraph:
            for vert1 in vertex1.edges:
                if vert1.mapping is not None:
                    if vert1.mapping not in vertex2.edges: return False

        # All of our tests have been passed, so the two vertices are a feasible pair
        return True

    ################################################################################

    def sort_label(self, dict p_feasible, dict init_map):
        return lambda vertex2: 0 if init_map is not None and vertex2 in init_map.values() else p_feasible

    cpdef preprocess(self, Graph graph1, Graph graph2, dict init_map, bint subgraph):
        """
        Sorts the vertices of the smaller graph into the order they will be searched in,
        then uses this ordering to pre-determine the terminal sets of the smaller graph.
        """

        cdef Vertex vertex1, vertex2, neighbor
        cdef dict label_counts = {}
        cdef dict degree_counts = {}

        cdef int graph1_size

        cdef int degree, max_degree

        cdef float p_label, p_degree, prob

        cdef dict p_feasible = {}

        cdef list exploration_sequence

        cdef int level, index

        cdef list terminal_levels = []  #at each call depth, the list of nodes that enter the terminal set
        cdef list core_ordering = []  #the order the nodes from graph2 will enter the mapping

        cdef dict parent_map = {}  # Map from each node to its parent

        graph1_size = len(graph1.vertices)

        # Count the number of times each label and degree appears in the larger graph
        for vertex1 in graph1.vertices:
            label_counts[vertex1.symbol] = label_counts.get(vertex1.symbol, 0) + 1
            degree_counts[len(vertex1.edges)] = degree_counts.get(len(vertex1.edges), 0) + 1

        # Find the probability of a node with a compatible label and degree
        # appearing in graph1 for each node in graph2
        if subgraph:
            # Find the maximum degree so we know what degrees we need to consider
            # when calculating probabilities on the vertices of the subgraph
            max_degree = 0

            for degree in degree_counts.keys():
                if degree > max_degree:
                    max_degree = degree

        for vertex2 in graph2.vertices:
            p_label = label_counts.get(vertex2.symbol, 0) / graph1_size

            if subgraph:
                p_degree = 0
                for degree in range(len(vertex2.edges), max_degree + 1):
                    p_degree += degree_counts.get(degree, 0) / graph1_size
            else:
                p_degree = degree_counts.get(vertex2.label, 0) / graph1_size

            p_feasible[vertex2] = p_degree * p_label

        #Create a sorted list of graph2's vertices from low -> high probability
        #Nodes within the initial mapping are moved to the front of the sequence
        exploration_sequence = sorted(graph2.vertices, key=self.sort_label(p_feasible, init_map))

        #Using the exploration sequence, pre-process the terminal sets of graph2 for each level of the matching
        level = 0
        for vertex2 in exploration_sequence:
            #Initialize the level to the previous one
            if level == 0:
                terminal_levels.append([])
            else:
                terminal_levels.append(terminal_levels[level - 1])

            #Add neighbors of the node to the terminal set
            for neighbor in vertex2.edges:
                terminal_levels[level].append(neighbor)

            #Remove nodes in the mapping from the terminal set
            for index in range(level + 1):
                if exploration_sequence[index] in terminal_levels[level]:
                    terminal_levels[level].remove(exploration_sequence[index])

            level += 1

        #Store the results of the pre-processing
        self.parents = parent_map
        self.exploration_sequence = exploration_sequence
        self.terminal_levels = terminal_levels

        print(exploration_sequence)
        print(terminal_levels)

        return

    #This was the isomorphism method assuming VF3 inherited VF2 for code reuse
    '''
    cdef isomorphism(self, Graph graph1, Graph graph2, dict initialMapping, bint subgraph, bint findAll):
        """
        Evaluate the isomorphism relationship between graphs `graph1` and
        `graph2` with optional initial mapping `initialMapping`. If `subgraph`
        is ``True``, `graph2` is treated as a possible subgraph of `graph1`.
        If `findAll` is ``True``, all isomorphisms are found; otherwise only
        the first is found.
        """

        #Peform pre-processing
        self.preprocess(graph1, graph2, initialMapping, subgraph)

        #The rest of the setup before calling match is the same
        #super.isomorphism will then call VF3's new match function
        super(VF3, self).isomorphism(self, graph1, graph2, initialMapping, subgraph, findAll)
    '''

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
            if vertex2.symbol is not vertex1.symbol:
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

    def debugPrinter(self):
        map = {}
        unmap1 = []
        unmap2 = []
        for vertex1 in self.graph1.vertices:
            if vertex1.mapping is not None:
                map[vertex1] = vertex1.mapping
            else:
                unmap1.append(vertex1)
        for vertex2 in self.graph2.vertices:
            if vertex2.mapping is None:
                unmap2.append(vertex2)

        print(map)
        print(unmap1)
        print(unmap2)