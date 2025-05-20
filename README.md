# Singular-simplicial-quiver-homologies


# Table of Contents

1. [Overview](#Overview)
3. [Installation](#Installation)  
4. [Usage](#Usage)
      1. [Example](#Example)
      3. [Reference manual](#Reference-manual)

# Overview

Extending singular homology to quivers in different categories, several homology theories have been constructed. We make use of isomorphisums between the quiver homologies and the homologies of certain spaces to created effective algorithms for computation.  Moreover, these constructions are natural and can be applied with persistent homology.


# Installation

Please note that the due to the use of the *Persistent Homology Algorithm Toolkit (phat)* for persistent homology computations with cell complexes, an earlier version of python and of some packages is required.

Required python packages: numpy (version 1.15.2), graphviz (version 0.14.2), phat (version 1.5.0).

Place a copy of "SingularQuiverHomologies.py" in your python path, nothing else is required other than Python 3 (version 3.5.6).


# Usage

Load the *SingularQuiverHomologies* library.

```python
import SingularQuiverHomologies as SQH
```


## Example

This example demonstrates the computation of all three singular simplicial homologies on a filtered quiver consisting of the a digraph suspension of a single double edge with an additional loop added at one of the vertices of the double edge. In particular, the flag complex is unaffected by the prescience of the loop and its homology coincides with that of a 2-sphere. Moreover, the reduced directed flag complex collapses the double edge. Hence is contactable at all filtration values. Finally, the reduced directed flag complex collapses the double edge once the loop enters the filtration. Therefore, has the homology of a 2-sphere until filtration value 0.25, after which it is contactable.

```python
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
```

The expected outputs of the above code are as follows.

<p align="center" name="Euler2">
      <img src="ConeQuiver.png" alt="alt text" width="100%" height="75%">
</p>

<p align="center" name="Euler2">
      <img src="ConeQuiverOutput.png" alt="alt text" width="100%" height="75%">
</p>


## Reference-manual



