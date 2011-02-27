"""
Roughly working toward Knuth's suggested solution for Ex 55.

A bunch of routines to generate ordered but NOT reduced
binary decision diagrams for the connectedness function
of a graph.

The idea is to run this and then feed the output into
Knuth's algorithm R in order to reduce it to a good old
friendly ordered reduced BDD.

Works okay for grids up to 8 by 8, but 9 by 9 starts to
get pretty ugly.
"""

import sys
import heapq

# define an ordering for the vertices by BFS from some root
def order_vertices(vertices, edges, root = None):
    closed = set()
    vertices = set(vertices)
    if root is None:
        root = vertices.pop()
    open = [(0, root)]
    ordering = []
    while open:
        d, vertex = heapq.heappop(open)
        if vertex not in closed:
            closed.add(vertex)
            ordering.append(vertex)
        for adj_vertex in edges[vertex]:
            if adj_vertex not in closed:
                heapq.heappush(open, (d + 1, adj_vertex))
    return ordering

# define an ordering for the edges based on the vertex ordering
def order_edges(vertices, edges, vertex_ordering):
    inverse_ordering = {}
    for (i, v) in enumerate(vertex_ordering):
        inverse_ordering[v] = i
    edge_ordering = []
    for u_i, u in enumerate(vertex_ordering):
        for v_i in sorted(inverse_ordering[v] for v in edges[u]):
            if v_i <= u_i:
                continue
            edge_ordering.append((u_i, v_i))
    return edge_ordering

# define the frontier sets based on the edge ordering
def make_frontiers(vertex_ordering, edge_ordering):
    n_vertices = len(vertex_ordering)
    frontiers = []
    # XXX TODO this might be wrong!
    for (u, v) in edge_ordering:
        frontiers.append(set(range(u, v + 1)))
    return frontiers

def make_counter():
    def gen_integers():
        i = 0
        while True:
            yield i
            i += 1
    integers = gen_integers().__iter__()
    def counter():
        return integers.next()
    return counter

def make_connectedness_tree(vertex_order, edge_order, frontiers, verbose = False):
    n_vertices = len(vertex_order)
    n_edges = len(edge_order)

    # we maintain a family of partitions of the vertices V in the graph
    # ie, a dict of lists of subsets of V, keyed by a unique index.
    # We'll increase the index as we generate the partitions
    # which we'll increase as we generate the partitions. this can be
    # processed later to form the bead keys for our BDD
    true_sink_index = -2
    false_sink_index = -1

    make_partition_index = make_counter()

    partitions = {make_partition_index():[set([0])]}
    beads = {}

    def make_bead(index, variable, low_index, high_index):
        beads[index] = (variable, low_index, high_index)

    make_bead(true_sink_index, n_edges, true_sink_index, true_sink_index)
    make_bead(false_sink_index, n_edges, false_sink_index, false_sink_index)

    def cached_partition(cache, partition, next_partitions, next_frontier_low):
        """
        add partition to cache, return index, frozen_partition
        """
        frozen_partition = tuple([
            frozenset(subset) for subset in partition
        ])
        if frozen_partition in cache:
            index = cache[frozen_partition]
            next_partitions[index] = frozen_partition
            return index
        elif len(partition) == 1 and len(partition[0]) == n_vertices:
            return true_sink_index
        elif next_frontier_low == n_vertices or any(max(s) < next_frontier_low for s in partition):
            return false_sink_index
        else:
            index = make_partition_index()
            cache[frozen_partition] = index
            next_partitions[index] = frozen_partition
            return index

    for depth, (edge, frontier) in enumerate(zip(edge_order, frontiers)):
        if verbose:
            print 'depth %d: beads %d, partitions %d' % (
                depth,
                len(beads),
                len(partitions)
            )
        if depth + 1 < n_edges:
            next_frontier_low = edge_order[depth + 1][0]
        else:
            next_frontier_low = n_vertices
        # cache partitions generated for each depth
        # this avoids a heap of duplication
        partition_cache = {}
        next_partitions = {}
        # branch on decision to include this edge
        for index, partition in partitions.iteritems():
            # Low subtree: don't include the edge
            low_partition = list(partition)
            if not any(edge[1] in subset for subset in partition):
                # add destination vertex as a singleton
                # subset if we haven't seen it yet
                low_partition.append(set([edge[1]]))

            low_index = cached_partition(
                partition_cache,
                low_partition,
                next_partitions,
                next_frontier_low,
            )

            # High subtree: include the edge
            high_partition = []
            to_connect = [set(edge)]
            for subset in partition:
                if edge[0] in subset or edge[1] in subset:
                    to_connect.append(subset)
                else:
                    high_partition.append(subset)
            newly_connected = set()
            for subset in to_connect:
                newly_connected = newly_connected.union(subset)
            high_partition.append(newly_connected)

            high_index = cached_partition(
                partition_cache,
                high_partition,
                next_partitions,
                next_frontier_low,
            )

            # create the bead for this split
            make_bead(index, depth, low_index, high_index)

        partitions = next_partitions

    """
    # evaluate each final partition for connectedness and set sink nodes
    for index, partition in partitions.iteritems():
        if len(partition) == 1 and len(partition[0]) == n_vertices:
            make_bead(index, depth + 1, true_sink_index, true_sink_index)
        else:
            make_bead(index, depth + 1, false_sink_index, false_sink_index)
    """
    # post-process - need to fix up the bead indices to agree with the usual
    # BDD indexing convention of root being highest, ..., 1 and 0 being the
    # true and false sinks
    s = len(beads)
    r = {
        -2 : 1,
        -1 : 0,
    }
    for i in xrange(s - 2):
        r[i] = s - 1 - i
    relabled_beads = {}
    def relabel(i):
        v, l, h = beads[i]
        relabled_beads[r[i]] = (v, r[l], r[h])
        del beads[i]
    relabel(-2)
    relabel(-1)
    for i in xrange(s - 2):
        relabel(i)
    return relabled_beads

def reduce_beads(beads):
    s = len(beads)
    redirect = {}
    cache = {}
    print 'reduce : making layers'

    layers = {}
    for (key, (v, l, r)) in beads.iteritems():
        if v not in layers:
            layers[v] = []
        layers[v].append((key, (v, l, r)))
    layer_variables = sorted(layers.keys(), reverse = True)

    def get_repr(key):
        while key in redirect:
            next_key = redirect[key]
            if next_key == key:
                break
            else:
                key = next_key
        return key

    print 'reduce : building redirects'
    # iterate over beads in decreasing variable order
    for variable in layer_variables:
        layer = layers[variable]
        for key, (v, l, r) in layer:
            l = get_repr(l)
            r = get_repr(r)
            if l == r:
                redirect[key] = l
            elif (v, l, r) in cache:
                redirect[key] = cache[(v, l, r)]
            else:
                cache[(v, l, r)] = key

    print 'reduce : walking'
    s = len(beads)
    root_key = s - 1
    reduced_beads = {}

    def walk(key):
        key = get_repr(key)
        if key in reduced_beads:
            return
        (v, l, r) = beads[key]
        l = get_repr(l)
        r = get_repr(r)
        reduced_beads[key] = (v, l, r)
        if l != key:
            walk(l)
        if r != key:
            walk(r)

    walk(root_key)
    print 'reduce : finished, returning'
    #XXX TODO fix keying? sigh.
    return rekey_monotone(reduced_beads)

def rekey_monotone(beads):
    r = {}
    for i, j in enumerate(sorted(beads.keys())):
        r[j] = i
    relabled_beads = {}
    def relabel(key):
        v, l, h = beads[key]
        relabled_beads[r[key]] = (v, r[l], r[h])
    for key in beads:
        relabel(key)
    return relabled_beads
