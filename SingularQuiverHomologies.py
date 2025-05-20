import numpy as np
from graphviz import Digraph
from copy import deepcopy
import phat
from itertools import combinations, chain
from multiprocessing import Pool


###General quiver functions###


#get vertices of quiver from edges
def Edge_list_to_quiver(edges):
    unfiltered_vertives = set([v for (u,v,f) in edges] + [v for (v,u,f) in edges])
    vertices = []
    for v in unfiltered_vertives:
        filtration_vlaues = [f for [u,w,f] in edges if (u == v or w == v)] 
        vertices.append([v, min(filtration_vlaues)])
    return [vertices, edges]


#removes all self loops form a list of edges
def remove_loops(edges):
    return [[u,w,f] for [u,w,f] in edges if u != w]


#records the vertices that have a loop and their corresponding minimal loop filtration values
def loop_vertex_filt_values(edges):
    loop_vertices = list(set([u for [u,w,f] in edges if u == w]))
    vertex_loop_min_filt_val = [min([f for [a,b,f] in edges if a == u and b == u]) for u in loop_vertices]
    return loop_vertices, vertex_loop_min_filt_val


#removes all multiple edges form a list of edges       
def remove_multiple_edges(edges):
    new_edges =  [[u,w,f] for [u,w,f] in edges if f == min([g for [a,b,g] in edges if (a == u and b == w)])]
    new_non_duplicated_edges = []
    [new_non_duplicated_edges.append(edge) for edge in new_edges if edge not in new_non_duplicated_edges]
    return new_edges


#removes loops and the multiple edges form a quiver whose filtration value are above a loop filtration value at either end
def remove_multiple_edges_above_loop_compute_edge_triples(edges, loop_vert, loop_vert_vals):
    new_edges = []
    loop_simplices = []
    non_loop_simplices = []
    for [u,v,f] in edges:
        #case when neither u nor v have a loop
        if u not in loop_vert and v not in loop_vert:
            loop_simplices.append(([u,v],[False,False]))
            non_loop_simplices.append(None)
            new_edges.append([u,v,f])
        #case when the edge is the first in the list to be added with a minimal filtration value
        elif [u, v, f] not in new_edges and f == min([g for [a,b,g] in edges if (a == u and b == v)]):
            new_edges.append([u,v,f])
            if u in loop_vert:
                if v in loop_vert:
                    loop_simplices.append(([u, v],[True, True]))
                    non_loop_simplices.append(-1)
                else:
                    loop_simplices.append(([u,v],[True, False]))
                    non_loop_simplices.append(0)
            else:
                loop_simplices.append(([u, v],[False, True]))
                non_loop_simplices.append(1)
        #otherwise the edge has a loop but is not added unless its filtration value is smaller than those of the incident loop(s)
        else:
            if u in loop_vert:
                if f < loop_vert_vals[loop_vert.index(u)]:
                    if v in loop_vert:
                        if f < loop_vert_vals[loop_vert.index(v)]:
                            new_edges.append([u,v,f])
                            loop_simplices.append(([u, v],[True, True]))
                            non_loop_simplices.append(-1)
                    else:
                        new_edges.append([u,v,f])
                        loop_simplices.append(([u, v],[True, False]))
                        non_loop_simplices.append(0)
            elif v in loop_vert:
                if f < loop_vert_vals[loop_vert.index(v)]:
                    new_edges.append([u,v,f])
                    loop_simplices.append(([u, v],[False, True]))
                    non_loop_simplices.append(1)
    return [new_edges, loop_simplices, non_loop_simplices]


#split quiver edges into those ending at a particular edge
def Split_edges_by_terminal_vertices(quiver):
    vertices = quiver[0]
    edges = quiver[1]
    vertex_terminal_edges = []
    for v in vertices:
        vertex_terminal_edges.append([[s,t,f] for [s,t,f] in edges if t == v[0]])#could remove found edges for grater efficency???
    return vertex_terminal_edges


#display function for filtered quiver
def Display_quiver(quiver, file_name = 'quiver.gv.pdf'):
    if file_name[-7:] != '.gv.pdf':
        file_name += '.gv.pdf'
    vertices = quiver[0]
    edges = quiver[1]
    quiver_plot = Digraph('Quiver',
                          filename = file_name ,
                          format="png")
    for vertex in vertices:
        quiver_plot.node(str(vertex[0]), label = str(vertex[0])+', '+str(vertex[1]))
    for edge in edges:
        quiver_plot.edge(str(edge[0]), str(edge[1]), label = str(edge[2]))
    quiver_plot.view()


###Directed flag complex functions###


#Extends k-partial (n+1)-dimensional simplices with maximal vertex v
#to k-partial (n+1)-dimensional simplices with maximal vertex v
def Extend_to_simplex(base_simplex, v_index_shift, n_simplices_at_v, previous_partial_simplices_at_v, k):
    #add base simplex to v-simplices and adjust irreverent loops up by one (or new if below)
    next_partial_simplices_at_v = []
    for partial_simplex in previous_partial_simplices_at_v:
        for sim_index in range(len(n_simplices_at_v)):
            possible_extension = True
            if n_simplices_at_v[sim_index][0] != base_simplex[k]:
                possible_extension = False
            else:
                for i in range(1, k + 1):
                    if n_simplices_at_v[sim_index][i] != n_simplices_at_v[partial_simplex[i] - v_index_shift][k]:
                        possible_extension = False
                        break
            if possible_extension:
                #sort filtration values
                filt_val = max(partial_simplex[-1], n_simplices_at_v[sim_index][-1])
                new_partial_simplex = deepcopy(partial_simplex)
                new_partial_simplex[-1] = sim_index + v_index_shift
                new_partial_simplex.append(filt_val)
                next_partial_simplices_at_v.append(new_partial_simplex)
    return next_partial_simplices_at_v


#Extends the n-skeleton of the directed flag complex of a quiver to the (n+1)-skeleton locally at a given top vertex v
def Extend_to_simplices_at_vertex(n, v, index_shifts, vertex_indexd_n_simplices):
    next_v_simplices = []
    if len(vertex_indexd_n_simplices[v]) != 0:
        num_vertices = len(vertex_indexd_n_simplices)
        for u in range(num_vertices):
            if u != v:
                for simplex_index in range(len(vertex_indexd_n_simplices[u])):
                    base_simplex = vertex_indexd_n_simplices[u][simplex_index]
                    new_v_simplices_from_simplex = [[simplex_index + index_shifts[u], base_simplex[-1]]]
                    for k in range(n):
                        new_v_simplices_from_simplex = Extend_to_simplex(base_simplex,
                                                                         index_shifts[v],
                                                                         vertex_indexd_n_simplices[v],
                                                                         new_v_simplices_from_simplex,
                                                                         k)
                    next_v_simplices += new_v_simplices_from_simplex
    return next_v_simplices


#Computes the directed flag complex of a quiver
def Directed_flag_complex(quiver, max_dim = 4):
    flag_complex = []
    if len(quiver) < 2:
        quiver = Edge_list_to_quiver(deepcopy(quiver))
    vertices = quiver[0]
    quiver[1] = remove_loops(quiver[1])
    v_indexed_simplices = Split_edges_by_terminal_vertices(quiver)
    edges = list(chain.from_iterable(v_indexed_simplices))#ensures the edge order is consistent
    flag_complex.append(vertices)
    flag_complex.append(edges)
    #extend the simplices
    num_vertices = len(vertices)
    start_sim_index = num_vertices
    for dim in range(2, max_dim + 1):
        new_v_indexed_simplices = []
        index_shifts = [start_sim_index]
        for v in range(num_vertices - 1):
            index_shifts.append(index_shifts[-1] + len(v_indexed_simplices[v]))
        for v in range(num_vertices):#parallel processing could be used here
            new_v_indexed_simplices.append(Extend_to_simplices_at_vertex(dim, v, index_shifts, v_indexed_simplices))
        #collect all simplices
        start_sim_index += len(flag_complex[-1])
        flag_complex.append(list(chain.from_iterable(new_v_indexed_simplices)))
        v_indexed_simplices = new_v_indexed_simplices
    return flag_complex


###Reduced flag complex functions###


#orders the vertices of abstract simplices
def Order_simplices(simplices, vertices):
    return [[v[0] for v in vertices if v[0] in s[:-1]] + [s[-1]] for s in simplices]


#with inputs simplices consisting of vertex lists (assumed to be order),
#outputs list of unique simplices by set of vertices
#with the minimal filtration value on simplices with the same vertices
def Reduce_vertes_ordered_simplices(simplices, all_simplices):
    unique_sim = []
    [unique_sim.append(s) for s in simplices if s[-1] == min([t[-1] for t in all_simplices if t[:-1] == s[:-1]]) and s not in unique_sim]
    return unique_sim


#converts simplices formed as lists of vertices to simplices containing a list of boundary indices of the simplices in the previous dimension
def Abstract_simplices_to_delta_set(ab_simplices, previous_ab_simplices, start_index):
    previous_ab_simplices = [s[:-1] for s in previous_ab_simplices]
    bnd_indexed_simplices = [[start_index+previous_ab_simplices.index(x) for x in [[s[j] for j in range(len(s)-1) if i!=j] for i in range(len(s)-1)]]+[s[-1]] for s in ab_simplices]
    return bnd_indexed_simplices


#converts simplices formed as lists of vertices to simplices containing a list of boundary indices of the simplices in the previous dimension
#and includes trivial case
def Delta_simplices(v_simplices, full_previous_simplices, start_index):
    if len(v_simplices) == 0:
        return []
    else:
        return Abstract_simplices_to_delta_set(v_simplices, full_previous_simplices, start_index)


#finds all the n+1 dimensional simplices that can be obtained by combining some n simplex and a particular vertex v
def Extend_simplex_by_vertex(simplices, v, v_maximal_edges):
    new_simplices = []
    for s in simplices:
        found = True
        edge_filtration_vals = []
        for i in range(len(s) - 1):
            if s[i] in v_maximal_edges[0]:
                edge_filtration_vals.append(v_maximal_edges[1][v_maximal_edges[0].index(s[i])])
            else:
                found = False
                break
        if found:
            new_simplices.append(s[:-1] + [v] + [max([s[-1]] + edge_filtration_vals)])
    return new_simplices


#Computes the reduced directed flag complex of a quiver
def Reduced_directed_flag_complex(quiver, max_dim = 4, processes = 1):
    redued_flag_complex = []
    quiver = deepcopy(quiver) 
    if len(quiver) < 2:
        quiver = Edge_list_to_quiver(quiver)
    quiver[1] = remove_loops(quiver[1])
    quiver[1] = remove_multiple_edges(quiver[1])
    vertices = quiver[0]
    v_indexed_edges = Split_edges_by_terminal_vertices(quiver)
    v_indexed_simplices = deepcopy(v_indexed_edges)
    v_indexed_edges = [([u[0] for u in v_indexed_simplices[v]], [u[-1] for u in v_indexed_simplices[v]]) for v in range(len(vertices))]
    edges = list(chain.from_iterable(v_indexed_simplices))#ensures the edge ordering is consistent
    redued_flag_complex.append(vertices)
    ordered_edges = Order_simplices(edges, vertices)
    redued_flag_complex.append(Reduce_vertes_ordered_simplices(ordered_edges, ordered_edges))
    full_simplices = redued_flag_complex[-1]
    #extend the simplices by vertices inductively for each dimension
    start_sim_index = len(vertices)
    for dim in range(2, max_dim + 1):
        new_v_indexed_simplices = []
        for v in range(len(vertices)):#parallel processing could be used here
            new_v_indexed_simplices.append(Extend_simplex_by_vertex(full_simplices, v, v_indexed_edges[v]))
        #order the vertices of the simplices
        ord_simplices = []
        for v in range(len(vertices)):#parallel processing could be used here
            ord_simplices.append(Order_simplices(new_v_indexed_simplices[v], vertices))
        full_ord_simplices = list(chain.from_iterable(ord_simplices))
        #covert new simplices to indexed delta complex format for reduced flag complex by removed non-unique simplices
        reduced_ord_simplices = []
        for v in range(len(vertices)):#parallel processing could be used here
            reduced_ord_simplices.append(Reduce_vertes_ordered_simplices(ord_simplices[v], full_ord_simplices))
        delta_simplices = []
        for v in range(len(vertices)):#parallel processing could be used here
            delta_simplices.append(Delta_simplices(reduced_ord_simplices[v], full_simplices, start_sim_index))
        start_sim_index += len(redued_flag_complex[-1])
        redued_flag_complex.append(list(chain.from_iterable(delta_simplices)))
        full_simplices = list(chain.from_iterable(reduced_ord_simplices))
        v_indexed_simplices = new_v_indexed_simplices
    return redued_flag_complex


###Partially reduced flag complex functions###

#Extends the n-skeleton of the directed flag complex of a quiver to the (n+1)-skeleton along with triples indicating loop and non-loop sub-simplices
def Extend_directed_flag_complex_triples(vertex_indexd_n_triples, n, loop_vertex_index, start_sim_index = 0):
    num_vertices = len(vertex_indexd_n_triples)
    next_vertex_indexd_triples = []
    index_shifts = [start_sim_index]
    #obtain shift in indices needed between simplices with different maximal vertices in their given order
    for v in range(num_vertices - 1):
        index_shifts.append(index_shifts[-1] + len(vertex_indexd_n_triples[v][0]))
    for v in range(num_vertices):#parallel processing could be used here
        next_v_triples = [[], [], []]
        loop_vertex = loop_vertex_index[v]
        if len(vertex_indexd_n_triples[v][0]) != 0: 
            for u in range(num_vertices):
                if u != v:
                    for simplex_index in range(len(vertex_indexd_n_triples[u][0])):
                        base_simplex = vertex_indexd_n_triples[u][0][simplex_index]
                        base_loop_siplex = vertex_indexd_n_triples[u][1][simplex_index]
                        base_non_loop_simplex = vertex_indexd_n_triples[u][2][simplex_index]

                        #obtain the full extend simplices as if they where in the flag complex
                        new_v_simplices_from_simplex = [[simplex_index + index_shifts[u], base_simplex[-1]]]
                        for k in range(n):
                            new_v_simplices_from_simplex = Extend_to_simplex(base_simplex,
                                                                             index_shifts[v],
                                                                             vertex_indexd_n_triples[v][0],
                                                                             new_v_simplices_from_simplex,
                                                                             k)
                        #deduce the other part of the triple, loop/non-loop sub-simplex structure
                        if loop_vertex:
                            next_v_triples[0] += new_v_simplices_from_simplex
                            next_v_triples[1] += [(base_loop_siplex[0] + [v], base_loop_siplex[1] + [True]) for s in range(len(new_v_simplices_from_simplex))]
                            next_v_triples[2] += [n for s in range(len(new_v_simplices_from_simplex))]
                        else:
                            next_v_triples[0] += new_v_simplices_from_simplex
                            next_v_triples[1] += [(base_loop_siplex[0] + [v], base_loop_siplex[1] + [False]) for s in range(len(new_v_simplices_from_simplex))]
                            next_v_triples[2] += [base_non_loop_simplex for s in new_v_simplices_from_simplex]
        next_vertex_indexd_triples.append(next_v_triples)
    return next_vertex_indexd_triples, index_shifts


#find the collection shortest paths between pairings to equally sized sets of distinct vertices in a given tree
def tree_paths(vertices1, vertices2, tree_edges):
    path_edges = []
    for v in vertices1:
        unuesed_tree_edges = list(np.arange(len(tree_edges)))
        paths = []
        boundery_vertices = [v]
        v_unmatched = True
        #extend all current paths by adjacent edges
        while v_unmatched:
            new_paths = []
            next_boundery = []
            for i in range(len(unuesed_tree_edges)):
                if tree_edges[unuesed_tree_edges[i]][0] in boundery_vertices:
                    path = next((p for p in paths if tree_edges[unuesed_tree_edges[i]][0] in tree_edges[p[-1]]), None)
                    if tree_edges[unuesed_tree_edges[i]][1] in vertices2:
                        if len(paths) > 0:
                            path_edges += path
                        path_edges.append(unuesed_tree_edges[i])
                        vertices2.remove(tree_edges[unuesed_tree_edges[i]][1])
                        v_unmatched = False
                        break
                    next_boundery.append(tree_edges[unuesed_tree_edges[i]][1])
                    if len(paths) == 0:
                        new_paths.append([unuesed_tree_edges[i]])
                    else:
                        new_paths.append(path + [unuesed_tree_edges[i]])
                elif tree_edges[unuesed_tree_edges[i]][1] in boundery_vertices:
                    path = next((p for p in paths if tree_edges[unuesed_tree_edges[i]][1] in tree_edges[p[-1]]), None)
                    if tree_edges[unuesed_tree_edges[i]][0] in vertices2:
                        if len(paths) > 0:
                            path_edges += path
                        path_edges.append(unuesed_tree_edges[i])
                        vertices2.remove(tree_edges[unuesed_tree_edges[i]][0])
                        v_unmatched = False
                        break
                    next_boundery.append(tree_edges[unuesed_tree_edges[i]][0])
                    if len(paths) == 0:
                        new_paths.append([unuesed_tree_edges[i]])
                    else:
                        new_paths.append(path + [unuesed_tree_edges[i]])
            paths = new_paths
            boundery_vertices = next_boundery
    #reduce path edges mod 2
    reduced_path_edges = []
    for edge in path_edges:
        if edge in reduced_path_edges:
            reduced_path_edges.remove(edge)
        else:
            reduced_path_edges.append(edge)
    return reduced_path_edges


#computes the minimal filtration values for extra cells between all top dimensional simplices on a set of vertices
def minimal_extra_cells(vertices, vertex_triples, previous_extra_cells, extra_indices, old_cell_shifts, cell_shifts, extra_shift, index_shift2, loop_verts = None, loop_vert_vals = None):
    num_vertices = len(vertices)
    extral_cells = []
    all_index_shifts = []
    all_vert_triples = [[], [], []]
    index_shift2 += extra_shift
    cell_shifts = [shift + extra_shift for shift in cell_shifts]
    tree_edges = [[c[0]-old_cell_shifts, c[1]-old_cell_shifts] for c in previous_extra_cells]
    for v in range(num_vertices):
        all_index_shifts += list(cell_shifts[vertices[v]] + np.arange(len(vertex_triples[v][0])))
        all_vert_triples[0] += vertex_triples[v][0]
        all_vert_triples[1] += vertex_triples[v][1]
        all_vert_triples[2] += vertex_triples[v][2]
    #determine all simplices on particular vertices
    all_vertex_simplices = [i for i in range(len(all_vert_triples[0])) if set(all_vert_triples[1][i][0]) == set(vertices)]
    #find the simplices that eventually merge in the complex hence need extra cells identifying them
    while len(all_vertex_simplices) > 1:
        current_simplices = [all_vertex_simplices[0]]
        del all_vertex_simplices[0]
        a = current_simplices[0]
        for i in range(len(all_vertex_simplices) - 1, -1, -1):
            b = all_vertex_simplices[i]
            if all_vert_triples[2][a] == all_vert_triples[2][b]:
                if set(np.array(all_vert_triples[1][a][0])[all_vert_triples[1][a][1]]) == set(np.array(all_vert_triples[1][b][0])[all_vert_triples[1][b][1]]):
                    current_simplices.append(b)
                    del all_vertex_simplices[i]
        #determine extra cells by minimal spanning tree of weights on all possible cells between vertex simplices
        #applies uses Prim’s Algorithm and compute simplex edge cell values on the fly
        if len(current_simplices) > 1:
            used_simplices = [current_simplices[0]]
            unused_simplices = current_simplices[1:]
            for i in range(len(current_simplices) - 1):
                min_filt = np.inf
                for j in used_simplices:
                    for k in unused_simplices:
                        other_cell_indices = []
                        other_cell_filt_vals = []
                        #check in any additional cells form the previous dimension are also faces on the new cell
                        for t in range(len(previous_extra_cells)):
                            vertices1 = list(set(all_vert_triples[0][j][:-1]) - set(all_vert_triples[0][k][:-1]))
                            vertices2 = list(set(all_vert_triples[0][k][:-1]) - set(all_vert_triples[0][j][:-1]))
                            other_cell_indices = tree_paths(vertices1, vertices2, tree_edges)
                            other_cell_filt_vals = [previous_extra_cells[t][-1] for l in other_cell_indices]
                            other_cell_indices = [index_shift2 + extra_indices[l] for l in other_cell_indices]
                        #in the dim 1 case we need to check the loop filtration values instead and keep the minimal among them
                        if type(loop_verts) != type(None) and type(loop_vert_vals) != type(None):
                            if all_vert_triples[0][j][0] in loop_verts:
                                loop_val_1 = loop_vert_vals[loop_verts.index(all_vert_triples[0][j][0])]
                                if all_vert_triples[0][j][1] in loop_verts:
                                    loop_val_2 = loop_vert_vals[loop_verts.index(all_vert_triples[0][j][1])]
                                    other_cell_filt_vals = [min([loop_val_1, loop_val_2])]
                                else:
                                    other_cell_filt_vals = [loop_val_1]
                            elif all_vert_triples[0][j][1] in loop_verts:
                                other_cell_filt_vals = [loop_vert_vals[loop_verts.index(all_vert_triples[0][j][1])]]
                        cell_filt = max([all_vert_triples[0][j][-1], all_vert_triples[0][k][-1]] + other_cell_filt_vals)
                        if cell_filt < min_filt:
                            next_sim = k
                            extra_cell = [all_index_shifts[j], all_index_shifts[k]] + other_cell_indices + [cell_filt]
                            min_filt = cell_filt
                used_simplices.append(next_sim)
                unused_simplices.remove(next_sim)
                extral_cells.append(extra_cell)
    return (vertices, extral_cells)


#for a specified vertex, evaluate stored face maps inductively to get the index on the maximal non-loop simplices in each dimension
def recover_non_loop_faces(v_simplices_n, v_faces_n, v_indexed_triples_n_1, shift_n_1, shift_n):
    num_vertices = len(v_indexed_triples_n_1)
    non_loop_faces_n_1 = []
    for u in range(num_vertices):
        non_loop_faces_n_1 += v_indexed_triples_n_1[u][2]
    new_v_faces_n = []
    for i in range(len(v_faces_n)):
        if type(v_faces_n[i]) == type(None):
            new_v_faces_n.append(shift_n + i)
        elif v_faces_n[i] == -1:
            new_v_faces_n.append(-1)
        else:
            new_v_faces_n.append(non_loop_faces_n_1[v_simplices_n[i][len(v_simplices_n[i]) - 2 - v_faces_n[i]] - shift_n_1])
    return new_v_faces_n


#replace None type with correctly indexed non-loop index when whole edge has no loops
def insert_missing_non_loop_edge_indices(v_sims, v_non_loop_sims, edge_shift):
    for i in range(len(v_non_loop_sims)):
        if type(v_non_loop_sims[i]) == type(None):
            v_non_loop_sims[i] = i + edge_shift
        elif v_non_loop_sims[i] >= 0:
            v_non_loop_sims[i] = v_sims[i][1 - v_non_loop_sims[i]]
    return v_non_loop_sims


#Computes the partially reduced directed flag complex of a quiver
def Partially_reduced_directed_flag_complex(quiver, max_dim = 4, processes = 1):

    quiver = deepcopy(quiver)

    #if only edges are supplied get vertices too
    if len(quiver) < 2:
        quiver = Edge_list_to_quiver(quiver)
    vertices = quiver[0]
    num_vertices = len(vertices)
        
    #get lists of vertex indices with loops along with their minimal loop filtration values
    loop_vertices, loop_vert_filt_vals = loop_vertex_filt_values(quiver[1])
    loop_vertex_index = [v in loop_vertices for [v,f] in vertices]
    
    #remove all loops from the quiver
    quiver[1] = remove_loops(quiver[1])
    
    #split edges by maximal vertices and edges
    v_indexed_edges = Split_edges_by_terminal_vertices(quiver)
    
    #remove all multiple edges with a loop filtration value lower than them
    v_indexed_simplex_triples = [[], []]
    for v in range(num_vertices):#parallel processing could be used here
        v_indexed_simplex_triples[1].append(remove_multiple_edges_above_loop_compute_edge_triples(v_indexed_edges[v], loop_vertices, loop_vert_filt_vals))#, start_index))
    edge_shift_indices = [len(vertices)]

    #initiate values of key variables
    partial_flag_complex = []
    partial_flag_complex.append(vertices)
    global_index_shifts = [0, len(vertices)]
    local_index_shifts = [[]]
    start_sim_index = len(vertices)
    
    #computes extend full complex of triples
    start_sim_index = num_vertices
    for dim in range(2, max_dim + 1):
        new_v_indexed_simplex_triples, dim_index_shifts = Extend_directed_flag_complex_triples(v_indexed_simplex_triples[-1], dim, loop_vertex_index, start_sim_index)
        start_sim_index += sum([len(v_indexed_simplex_triples[-1][v][0]) for v in range(num_vertices)])
        global_index_shifts.append(start_sim_index)
        local_index_shifts.append(dim_index_shifts)
        v_indexed_simplex_triples.append(new_v_indexed_simplex_triples)
    local_index_shifts.append([start_sim_index])
    for v in range(num_vertices - 1):
        local_index_shifts[-1].append(local_index_shifts[-1][-1] + len(new_v_indexed_simplex_triples[v][0]))
    
    #evaluate stored face maps inductively to get the index on the maximal non-loop simplices in each dimension
    for v in range(num_vertices):#parallel processing could be used here
        v_indexed_simplex_triples[1][v][2] = insert_missing_non_loop_edge_indices(v_indexed_simplex_triples[1][v][0],
                                                                                  v_indexed_simplex_triples[1][v][2],
                                                                                  local_index_shifts[1][v])
    for dim in range(2, max_dim + 1):
        for v in range(len(vertices)):#parallel processing could be used here
            v_indexed_simplex_triples[dim][v][2] = recover_non_loop_faces(v_indexed_simplex_triples[dim][v][0], 
                                                                          v_indexed_simplex_triples[dim][v][2],
                                                                          v_indexed_simplex_triples[dim - 1],
                                                                          global_index_shifts[dim - 1],
                                                                          local_index_shifts[dim][v])
    
    #computes full partial flag complex required for persistence computations form the triples complex by adding in additional cells when required
    extra_index_shift = 0
    old_extra_index_shift = 0
    num_last_cells = 0
    previous_extra_cells = []
    extra_indices = []
    all_vertices = list(np.arange(num_vertices))
    for dim in range(1, max_dim + 1):
        all_dim_simplices = [sim for v in range(num_vertices) for sim in v_indexed_simplex_triples[dim][v][0]]
        if dim > 1:
            partial_flag_complex.append([[s[i] + extra_index_shift for i in range(len(s)-1)] + [s[-1]] for s in all_dim_simplices] + extra_cells)
            old_extra_index_shift = extra_index_shift
            extra_index_shift += num_last_cells
            num_last_cells = len(extra_cells)
            dim_1_loop_vert = None,
            dim_1_loop_vert_vals = None
        else:
            partial_flag_complex.append([[s[i] for i in range(len(s)-1)] + [s[-1]] for s in all_dim_simplices])
            dim_1_loop_vert = loop_vertices
            dim_1_loop_vert_vals = loop_vert_filt_vals
        extra_cells = []
        if dim != max_dim:
            total_dim_simplices = sum([len(v_sims[0]) for v_sims in v_indexed_simplex_triples[dim]])
            for vertices in combinations(all_vertices, dim + 1):#parallel processing could be used here
                extra_cells.append(minimal_extra_cells(vertices,
                                                    [v_indexed_simplex_triples[dim][v] for v in vertices],
                                                    [c for (v, cells) in previous_extra_cells for c in cells if set(v).issubset(vertices)],
                                                    [index for c in range(len(extra_indices)) for index in extra_indices[c] if set(previous_extra_cells[c][0]).issubset(vertices)],
                                                    old_extra_index_shift,
                                                    local_index_shifts[dim],
                                                    extra_index_shift,
                                                    global_index_shifts[dim] + total_dim_simplices,
                                                    dim_1_loop_vert,
                                                    dim_1_loop_vert_vals))
            previous_extra_cells = deepcopy(extra_cells)
            extra_indices = []
            count = 0
            for cells in extra_cells:
                extra_indices.append(list(range(count, count + len(cells[1]))))
                count += len(cells[1])
            extra_cells = [cell for (v, cells) in extra_cells for cell in cells]
    return partial_flag_complex


###Persistent homology functions###


#concatenate simplices and sort filter values
def Filter_ordered_simplices(space):
    all_dimensions = list(chain.from_iterable([[d for i in range(len(space[d]))] for d in range(len(space))]))
    all_simplices = list(chain.from_iterable(space))
    sort_index = sorted(range(len(all_simplices)), key = lambda i:all_simplices[i][-1])
    ordered_simplicies = deepcopy([[sort_index.index(x) for x in all_simplices[sort_index[i]][:-1]] + [all_simplices[sort_index[i]][-1]] for i in range(len(sort_index))])
    dimensions = [all_dimensions[i] for i in sort_index]
    return ordered_simplicies, dimensions

#compute persistent homology pairs in each dimension
def Persistent_homology(space):
    filter_vlaues = []
    dimensions = []
    all_simplices, dimensions = Filter_ordered_simplices(space)
    #form filtered boundary matrix
    boundary_matrix = phat.boundary_matrix(representation = phat.representations.vector_vector)
    columns = []
    for s in range(len(all_simplices)):
        filter_vlaues.append(all_simplices[s][-1])
        if dimensions[s] == 0:
            columns.append((0, []))
        else:
            columns.append((dimensions[s], all_simplices[s][:-1]))
    boundary_matrix.columns = columns
    #compute homology
    pairs = boundary_matrix.compute_persistence_pairs_dualized(reduction=phat.reductions.chunk_reduction)
    #collect pairs by dimension and with true birth and death values
    pairs_by_dim = []
    inf_points_by_dim = []
    for i in range(len(space)):
        inf_points_by_dim.append([])
        pairs_by_dim.append([])
    inf_index = list(range(len(all_simplices)))
    for i in range(len(pairs)):
        birth = filter_vlaues[pairs[i][0]]
        inf_index.remove(pairs[i][0])
        death = filter_vlaues[pairs[i][1]]
        inf_index.remove(pairs[i][1])
        #remove pairs if they have the same birth and death filtration values
        if birth != death:
            pairs_by_dim[dimensions[pairs[i][0]]].append([birth, death])
    for i in inf_index:
        inf_points_by_dim[dimensions[i]].append(all_simplices[i][-1])
    return pairs_by_dim, inf_points_by_dim

