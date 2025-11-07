#!/bin/bash

## tc_grps for pulling simulation: indexing for groups
gmx make_ndx -f step7_production.gro -o index.ndx <<EOF
1 | r BQ1 | r BQ2 | r BQ3 | r BQ4 | r BQ5 | r BQ6 | r BQ7
name 24 Protein_BQ
r TIP3 | r SOD | r CLA
name 25 W_Ion
r POPC | r TIP3 | r SOD | r CLA
name 26 POPC_W_Ion
a 281  752  764  766  807  808  810  812  813  814  820  823  844  845  852 \
 853  854  855  862  863  864  866  867  868  869  873  922  927 1896 1897 \
1898 1900 1938 1939 1942 1943 1945 1946 1947 1948 2001 2006 2007 2008 2009 \
12677 12679 12680 12701 12703 12705 12706 12707 12710 12721 12722 12723 12769 12770 12772
name 27 3A_within_Protein
q
EOF

## create group of protein atoms within 3 Å (0.3 nm) of BQ2
gmx select -s step7_production.tpr -n index.ndx \
-select 'group "Protein" and within 0.3 of group "BQ1"' \
-on index_within3A_BQ2.ndx

### add this index_within3A_BQ2.ndx to index.ndx manually

## pulling simulation
gmx grompp -f md_pull.mdp -c step7_production.gro -r step7_production.gro \
-n index.ndx -p system.top -o pull.tpr

gmx mdrun -deffnm pull -nt 10 -gpu_id 0 -v

