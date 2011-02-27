"""
Plot uniformly sampled random connected subgraphs of n*n grid.
"""

import numpy

from connection import order_vertices, order_edges, make_frontiers, \
    make_connectedness_tree, reduce_beads

from bdd import bdd_count_solutions, bdd_generate_random_solution

def make_grid_graph(n):
    """
    makes graph (V, E) of n^2 grid
    """
    vertices = set()
    edges = {}
    for i in xrange(n):
        for j in xrange(n):
            vertices.add((i, j))

    for (i, j) in vertices:
        for (i2, j2) in ((i, j-1), (i, j+1), (i-1, j), (i+1, j)):
            if (i2, j2) in vertices:
                edges[(i, j)] = edges.get((i, j), []) + [(i2, j2)]
    return (vertices, edges)

def make_bmp_graph(bmp):
    assert len(bmp.shape) == 2
    n, m = bmp.shape
    vertices = set()
    edges = {}
    for i in xrange(n):
        for j in xrange(m):
            if bmp[i, j]:
                vertices.add((i, j))

    for (i, j) in vertices:
        for (i2, j2) in ((i, j-1), (i, j+1), (i-1, j), (i+1, j)):
            if (i2, j2) in vertices:
                edges[(i, j)] = edges.get((i, j), []) + [(i2, j2)]
    return (vertices, edges)

def gen_random_solutions(beads, how_many):
    bdd_beads = {
        's' : len(beads),
        'dag' : beads,
    }
    c = {}
    n_solns = bdd_count_solutions(bdd_beads, c)
    print 'number of solutions : %d' % n_solns
    print 'here are a few random ones:'
    for _ in xrange(how_many):
        yield bdd_generate_random_solution(
            bdd_beads,
            c,
            rand = numpy.random.rand
        )

def make_bdd(vertices, edges, root):
    vertex_order = order_vertices(vertices, edges, root)
    edge_order = order_edges(vertices, edges, vertex_order)
    frontiers = make_frontiers(vertex_order, edge_order)

    print 'begin horrific connectedness tree construction procedure'
    beads = make_connectedness_tree(
        vertex_order,
        edge_order,
        frontiers,
        verbose = True,
    )
    print '\noutput (unreduced) contains %d beads\n' % len(beads)
    beads = reduce_beads(beads)
    print '\noutput (reduced) contains %d beads\n' % len(beads)
    return beads, edge_order, vertex_order

def make_bmp(n, edge_order, vertex_order, soln):
    bmp = numpy.zeros((2 * n + 1, ) * 2, dtype = numpy.int)

    def carve_edge(bmp, edge_i):
        (u_i, v_i) = edge_order[edge_i]
        (u_x, u_y) = vertex_order[u_i]
        (v_x, v_y) = vertex_order[v_i]
        bmp[2 * u_x + 1, 2 * u_y + 1] = 1
        bmp[2 * v_x + 1, 2 * v_y + 1] = 1
        if u_x == v_x and u_y == v_y + 1:
            bmp[2 * u_x + 1, 2 * v_y + 2] = 1
        elif u_x == v_x and u_y + 1 == v_y:
            bmp[2 * u_x + 1, 2 * u_y + 2] = 1
        elif u_x == v_x + 1 and u_y == v_y:
            bmp[2 * v_x + 2, 2 * u_y + 1] = 1
        elif u_x + 1 == v_x and u_y == v_y:
            bmp[2 * u_x + 2, 2 * u_y + 1] = 1

    for i, edge in enumerate(soln):
        if edge:
            carve_edge(bmp, i)
    return bmp

def coarsen_bmp(bmp, coarse_factor):
    coarse_shape = tuple(coarse_factor * n for n in bmp.shape)
    coarse_bmp = numpy.zeros(coarse_shape)
    for i in xrange(coarse_factor):
        for j in xrange(coarse_factor):
            coarse_bmp[i::coarse_factor, j::coarse_factor] = bmp
    return coarse_bmp

def print_bmp(bmp):
    for line in bmp:
        print ''.join(['#' if x else ' ' for x in line])

def main():
    # trying anything above n = 5 may prove a bit foolish
    n = 2
    print 'making a %d by %d grid' % (n, n)
    vertices, edges = make_grid_graph(n)
    # using a corner vertex as root works much better
    # than a central one in terms of reducing the
    # size of the connectedness tree
    beads, edge_order, vertex_order = make_bdd(vertices, edges, root = (0, 0))

    subplot_width = 2 * n + 1
    subplot_height = 2 * n + 1

    (soln, ) = list(gen_random_solutions(beads, how_many = 1))

    bmp = make_bmp(n, edge_order, vertex_order, soln)
    print_bmp(bmp)
    coarse_bmp = coarsen_bmp(bmp, coarse_factor = 2)
    print_bmp(coarse_bmp)

    c_vertices, c_edges = make_bmp_graph(coarse_bmp)
    print 'n vertices: %d; n edges: %d' % (len(c_vertices), len(c_edges))
    c_root = min(c_vertices)
    beads, edge_order, vertex_order = make_bdd(c_vertices, c_edges, c_root)

    for soln in gen_random_solutions(beads, how_many = 25):
        bmp = make_bmp(coarse_bmp.shape[0], edge_order, vertex_order, soln)
        print_bmp(bmp)

if __name__ == '__main__':
    main()
