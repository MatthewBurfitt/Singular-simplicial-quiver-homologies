import SingularQuiverHomologies as SQH

#double edge cones with loop
edges = [[0,1,0.0],[1,0,0.1],[0,2,0.2],[0,0,0.25],[1,2,0.3],[0,3,0.4],[1,3,0.5]]

quiver = SQH.Edge_list_to_quiver(edges)
SQH.Display_quiver(quiver, file_name = 'MultiEdgeConeQuiver.gv.pdf')


#flag complex
flag_complex = SQH.Directed_flag_complex(quiver, max_dim = 5)

print('\n Flag complex of the cone has simplices:')
for i in range(len(flag_complex)):
        print(len(flag_complex[i]), i, 'simplices')

pairs_by_dim, inf_points_by_dim = SQH.Persistent_homology(flag_complex)

print('\n Flag complex of the cone has persistence:')
print(' ')
for i in range(len(pairs_by_dim)):
        if len(pairs_by_dim[i]) != 0:
                print('number of persistence pairs in dim', i, ' is ', len(pairs_by_dim[i]))
        print('... and 0 otherwise')
        print(' ')
        for i in range(len(inf_points_by_dim)):
                if len(inf_points_by_dim[i]) != 0:
                    print('number of infinite points in dim', i, ' is ', len(inf_points_by_dim[i]))
print('... and 0 otherwise')


#reduced flag complex
reduced_flag_complex = SQH.Reduced_directed_flag_complex(quiver, max_dim = 5)

print('\n Reduced Flag complex of the cone has simplices:')
for i in range(len(reduced_flag_complex)):
        print(len(reduced_flag_complex[i]), i, 'simplices')

reduced_pairs_by_dim, reduced_inf_points_by_dim = SQH.Persistent_homology(reduced_flag_complex)

print('\n Reduced flag complex of the cone has persistence:')
print(' ')
for i in range(len(reduced_pairs_by_dim)):
        if len(reduced_pairs_by_dim[i]) != 0:
                print('number of persistence pairs in dim', i, ' is ', len(reduced_pairs_by_dim[i]))
print('... and 0 otherwise')
print(' ')
for i in range(len(reduced_inf_points_by_dim)):
        if len(reduced_inf_points_by_dim[i]) != 0:
                print('number of infinite points in dim', i, ' is ', len(reduced_inf_points_by_dim[i]))
print('... and 0 otherwise')


#partially reduced flag complex
partial_flag_complex = SQH.Partially_reduced_directed_flag_complex(quiver, max_dim = 5)

print('\n Partially reduced Flag complex of the cone has simplices:')
for i in range(len(partial_flag_complex)):
        print(len(partial_flag_complex[i]), i, 'simplices')

partial_pairs_by_dim, partial_inf_points_by_dim = SQH.Persistent_homology(partial_flag_complex)

print('\n Partially reduced flag complex of the cone has persistence:')
print(' ')
for i in range(len(partial_pairs_by_dim)):
        if len(partial_pairs_by_dim[i]) != 0:
                print('number of persistence pairs in dim', i, ' is ', len(partial_pairs_by_dim[i]))
print('... and 0 otherwise')
print(' ')
for i in range(len(partial_inf_points_by_dim)):
        if len(partial_inf_points_by_dim[i]) != 0:
                print('number of infinite points in dim', i, ' is ', len(partial_inf_points_by_dim[i]))
print('... and 0 otherwise')

