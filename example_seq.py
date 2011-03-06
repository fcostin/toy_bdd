"""
Plot uniformly sampled random connected subgraphs of n*n grid.
"""

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


def count_em(n):
    vertices, edges = make_grid(n)
    corner_root = (0, ) * 2
    vertex_order = order_vertices(vertices, edges, root = corner_root)
    edge_order = order_edges(vertices, edges, vertex_order)
    frontiers = make_frontiers(vertex_order, edge_order)
    beads = make_connectedness_tree(
        vertex_order,
        edge_order,
        frontiers,
    )
    beads = reduce_beads(beads, verbose = False)
    count = bdd_count_solutions({'s' : len(beads), 'dag' : beads})
    print count * 2

def main():
    for n in xrange(2, 5 + 1):
        count_em(n)

if __name__ == '__main__':
    main()
