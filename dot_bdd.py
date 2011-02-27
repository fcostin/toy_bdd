"""
Export a BDD graph to file for plotting with graphviz' dot

Beads ordered according to branch variable.

For example, after exporting your beads to 'bdd.gv' you
can run the following to layout and plot the graph:

    dot -v -Tpng -O bdd.gv
    eog bdd.gv.png
"""

def export_dot_graph(beads, file_name):
    out_file = open(file_name, 'w')
    write_line = lambda s : out_file.write(s + '\n')
    write_line('digraph tree {')
    write_line('\tgraph []')
    layers = {}
    for (key, (v, l, r)) in beads.iteritems():
        if v not in layers:
            layers[v] = []
        layers[v].append((key, (v, l, r)))
    for value in sorted(layers):
        layer_beads = layers[value]
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
                # n.b. html entity (decimal) for unicode uptack (upside down T)
                label = 'T' if index else '&#8869;'
                write_line('\t\t"%d" [label="%s", shape = box];' % (index, label))
            else:
                write_line('\t\t"%d" [label="%d", shape = circle];' % (index, value))
                write_line('\t\t"%d" -> "%d" [style=dashed];' % (index, left_index))
                write_line('\t\t"%d" -> "%d" [style=solid];' % (index, right_index))

        write_line('\t}')
    write_line('}')
    out_file.close()


