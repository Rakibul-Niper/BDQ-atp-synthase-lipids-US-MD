#!/bin/bash


gmx make_ndx -f us-windows/window0.gro -o index.ndx << EOF
1 | r BQ1 | r BQ2 | r BQ3 | r BQ4 | r BQ5 |r BQ6 | r BQ7 | r POPC | r TIP3 | r SOD | r CLA
name 24 Protein-BQ-POPC-SOD-CLA
1 | r BQ1 | r BQ2 | r BQ3 | r BQ4 | r BQ5 |r BQ6 | r BQ7 
name 25 Protein_BQ
r TIP3 | r SOD | r CLA
name 26 W_Ion
r POPC | r TIP3 | r SOD | r CLA
name 27 POPC_W_Ion
a 281  752  764  766  807  808  810  812  813  814  820  823  844  845  852 \
 853  854  855  862  863  864  866  867  868  869  873  922  927 1896 1897 \
1898 1900 1938 1939 1942 1943 1945 1946 1947 1948 2001 2006 2007 2008 2009 \
12677 12679 12680 12701 12703 12705 12706 12707 12710 12721 12722 12723 12769 12770 12772
name 28 3A_within_Protein
q
EOF


mkdir -p umbrella-md      
cd umbrella-md           
count=0
targets=( $(seq 0.50 0.10 1.50) $(seq 1.60 0.10 3.00) )

echo "📂 Starting Umbrella Sampling MD Runs..."
echo "Total windows: ${#targets[@]}"
echo "====================================="

for dist in "${targets[@]}"; do
    printf "\n🚀 Running window %02d at %.2f nm\n" "$count" "$dist"
    sed "s/XX/${dist}/g" ../inputs/us.template > ../inputs/md_umbrella${count}.mdp
    gmx grompp -f ../inputs/md_umbrella${count}.mdp \
               -c ../us-windows/window${count}.gro -r ../us-windows/window${count}.gro \
               -n ../index.ndx -p ../system.top \
               -o md_umbrella${count}.tpr -maxwarn 1
    gmx mdrun -v -nt 10 -gpu_id 0 -deffnm md_umbrella${count}
    printf "✅ Completed window %02d at %.2f nm\n" "$count" "$dist"
    ((count++))
done

echo "====================================="
echo "🎉 Umbrella sampling MD finished for $count windows!"


