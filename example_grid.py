"""
Plot uniformly sampled random connected subgraphs of n*n grid.
"""

import pylab
import numpy

from connection import order_vertices, order_edges, make_frontiers, \
    make_connectedness_tree, reduce_beads

from bdd import bdd_count_solutions, bdd_generate_random_solution

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

def main():
    # trying anything above n = 5 may prove a bit foolish
    n = 5
    print 'making a %d by %d grid' % (n, n)
    vertices, edges = make_grid(n)
    # experiment: trying to fix roots
    central_root = (n/2, ) * 2 # this seems to work poorly
    corner_root = (0, ) * 2
    vertex_order = order_vertices(vertices, edges, root = corner_root)
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

    plots_across = 5
    plots_down = 5
    n_plots = plots_across * plots_down

    subplot_margin = -1
    subplot_width = 2 * n + 1
    subplot_height = 2 * n + 1
    plot_width = (
        plots_across * subplot_width + subplot_margin * (plots_across - 1)
    )
    plot_height = (
        plots_down * subplot_height + subplot_margin * (plots_down - 1)
    )

    figure_bmp = numpy.ones((plot_width, plot_height), dtype = numpy.int)

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

    for plot, soln in enumerate(gen_random_solutions(beads, how_many = n_plots)):
        bmp = numpy.zeros((subplot_width, subplot_height), dtype = numpy.int)
        for i, include_edge in enumerate(soln):
            if include_edge:
                carve_edge(bmp, i)

        plot_x = plot / plots_down
        plot_y = plot % plots_down
        plot_xx = plot_x * (subplot_width + subplot_margin)
        plot_yy = plot_y * (subplot_height + subplot_margin)
        x_slice = slice(plot_xx, plot_xx + subplot_width)
        y_slice = slice(plot_yy, plot_yy + subplot_height)
        figure_bmp[x_slice, y_slice] = 4 * (1 - bmp)
        # print soln

    pylab.figure()
    pylab.imshow(
        figure_bmp.T,
        interpolation = 'nearest',
        cmap = pylab.cm.Greys,
    )
    pylab.xticks([])
    pylab.yticks([])
    pylab.imsave(
        'pretties.png',
        figure_bmp.T,
        cmap = pylab.cm.Greys,
    )
    pylab.show()

if __name__ == '__main__':
    main()
