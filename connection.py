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
from itertools import groupby

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
        # print 'making frontier; u, v = %d, %d' % (u, v)
        frontiers.append(set(range(u, v + 1)))
    return frontiers

def make_grid(n):
    """
    makes graph (V, E) of n by n grid.
    """
    vertices = set()
    edges = {}
    for i in xrange(n):
        for j in xrange(n):
            vertex = (i, j)
            vertices.add(vertex)
            for (i2, j2) in ((i, j-1), (i, j+1), (i-1, j), (i+1, j)):
                if i2 < 0 or i2 == n or j2 < 0 or j2 == n:
                    continue
                edges[vertex] = edges.get(vertex, []) + [(i2, j2)]
    return (vertices, edges)

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

def test():
    # trying anything above n = 5 may prove a bit foolish
    n = 4
    print 'making a %d by %d grid' % (n, n)
    vertices, edges = make_grid(n)
    # experiment: trying to fix roots
    # central_root = (n/2, ) * 2 # this seems to work poorly
    corner_root = (0, ) * 2
    vertex_order = order_vertices(vertices, edges, root = corner_root)
    edge_order = order_edges(vertices, edges, vertex_order)
    frontiers = make_frontiers(vertex_order, edge_order)
    for depth, (edge, frontier) in enumerate(zip(edge_order, frontiers)):
        print 'depth %d edge %s frontier %s' % (
            depth,
            str(edge),
            str(frontier),
        )

    print 'begin horrific connectedness tree construction procedure'
    beads = make_connectedness_tree(
        vertex_order,
        edge_order,
        frontiers,
        verbose = True,
    )
    print ''
    print 'output (unreduced) contains %d beads' % len(beads)
    print ''

    beads = reduce_beads(beads)
    print ''
    print 'output (reduced) contains %d beads' % len(beads)
    print ''
    if len(beads) <= 1000:
        dump_graph(beads, 'tree.gv')
        sys.exit(0)
    else:
        print 'refusing to dump graph as n too large'
        sys.exit(0)

def dump_graph(beads, file_name):
    out_file = open(file_name, 'w')
    write_line = lambda s : out_file.write(s + '\n')
    write_line('digraph tree {')
    write_line('\tgraph []')
    layers = {}
    for value, layer_iter in groupby(beads.iteritems(), key = lambda (k, (v, l, r)) : v):
        layers[value] = list(layer_iter)
    for value in sorted(layers):
        layer_beads = list(layers[value])
        write_line('\t{')
        write_line('\t\trank = same;')
        for index, bead in layer_beads:
            (value, left_index, right_index) = bead
            write_line('\t\t"%d" [label="%d"];' % (index, value))
        write_line('\t}')
        write_line('\t{')
        for index, bead in layer_beads:
            (value, left_index, right_index) = bead
            if left_index == index and right_index == index:
                is_sink = True
            elif left_index != index and right_index != index:
                is_sink = False
            else:
                raise ValueError('bad bead : %s' % str(bead))
            if is_sink:
                label = 'T' if index else '&#8869;' # html entity (decimal) for unicode uptack (upside down T)
                write_line('\t\t"%d" [label="%s", shape = box];' % (index, label))
            else:
                write_line('\t\t"%d" [label="%d", shape = circle];' % (index, value))
                write_line('\t\t"%d" -> "%d" [style=dashed];' % (index, left_index))
                write_line('\t\t"%d" -> "%d" [style=solid];' % (index, right_index))

        write_line('\t}')
    write_line('}')
    out_file.close()

def reduce_beads(beads):
    s = len(beads)
    redirect = {}
    cache = {}
    print 'reduce : making layers'
    layers = [(v, list(i)) for (v, i) in groupby(beads.iteritems(), key = lambda (k, (v, l, r)) : v)]
    layers = sorted(layers, reverse = True)

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
    for variable, layer in layers:
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
        (v, l, r) = beads[key]
        l = get_repr(l)
        r = get_repr(r)
        reduced_beads[key] = (v, l, r)
        if l != key:
            walk(l)
        if r != key:
            walk(r)
    walk(root_key)
    #XXX TODO  fix numbering? sigh.
    return reduced_beads

if __name__ == '__main__':
    test()
