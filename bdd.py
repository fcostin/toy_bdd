"""
toy BDD implementation based on Knuth's AoCP, Vol 4 Fascicle 1

BDD consists of a DAG of beads
each bead is identified with a key k and the value (k_v, k_l, k_h)
here, k_v is the index of the variable to test, and we jump to the key k_l next
if x_v is 0, and jump to k_h otherwise.

A few monotonicity requirements must be satisfied:

    1. k_l, k_h < k for all keys k.
    this ensures that we have a DAG.

    2. k_l_v, k_h_v > k_v for all keys k.
    this ensures that the variables in the
    BDD are tested in order.

The special beads for the sinks (True and False) break both of the above rules.
First, k_l = k_h = k for both sinks. This breaks the first condition, which implies
that the second condition cannot hold either. For the second condition to hold
for all non-sink beads in the BDD we need to assign keys k_true and k_false for
the sink nodes that are strictly less than all other keys in the structure. I
guess we can do this by convention, by always ensuring the True sink has key 1,
and the False sink has key 0.
"""

BDD_X = dict(
    n_args = 4,
    s = 9,
    dag = {8:(1,7,6), 7:(2,5,4), 6:(2,0,1), 5:(3,1,0), 4:(3,3,2), 3:(4,1,0), 2:(4,0,1), 1:(5,1,1), 0:(5,0,0)},
)

BDD_Y = dict(
    n_args = 4,
    s = 9,
    dag = {8:(1,7,2), 7:(2,4,6), 6:(3,3,5), 5:(4,0,1), 4:(3,1,0), 3:(4,1,0), 2:(2,0,1), 1:(5,1,1), 0:(5,0,0)},
)

def bdd_root(bdd):
    # nb this implies convention that s=2 of a 1 bead BDD containing only True sink (Knuth p75)
    return bdd['s'] - 1

"""
my soln to Ex. 10
"""
def bdd_equality(x, y):
    x_dag = x['dag']
    y_dag = y['dag']

    # do a DFS over the pair of dags checking each pair of nodes are equal
    # (either both the same sink, or both branching on the same variable)

    def compare(x_key, y_key):
        if x_dag[x_key][0] != y_dag[y_key][0]:
            return False
        elif (x_key < 2) or (y_key < 2):
            return x_key == y_key
        else:
            return compare(x_dag[x_key][1], y_dag[y_key][1]) and compare(x_dag[x_key][2], y_dag[y_key][2])

    return compare(bdd_root(x), bdd_root(y))

"""
algorithm C. p75, Knuth
"""
def bdd_count_solutions(bdd, c = None):
    """
    compute number of solutions of given bdd. If optional argument c is given,
    it should be a map, which will then be used to tabulate the counts of the
    sub graphs for each bead in the bdd.
    """
    s = bdd['s']
    dag = bdd['dag']
    if c is None:
        c = {}
    c[0] = 0
    c[1] = 1
    for k in xrange(2, s):
        v_k, l, h = dag[k]
        v_l = dag[l][0]
        v_h = dag[h][0]
        # nb the exponential weights are required because we can jump from testing
        # a variable x_i to a variable x_j with j > i + 1. this skips the variables
        # between i and j, so we mulitply the count by 2 for each one we skip..
        c[k] = (2 ** (v_l - v_k - 1)) * c[l] + (2 ** (v_h - v_k - 1)) * c[h]

    # first variable tested, ie by the root bead, may not be the first variable
    # so as above we multiply by 2 for each skipped variable
    v_root = dag[s - 1][0]
    return 2 ** (v_root - 1) * c[s - 1]

def main():
    print 'are BDDs x and y equal? %s' % bdd_equality(BDD_X, BDD_Y)

    print 'how many solutions does x have? %d' % bdd_count_solutions(BDD_X)
    print 'how many solutions does y have? %d' % bdd_count_solutions(BDD_Y)

if __name__ == '__main__':
    main()
