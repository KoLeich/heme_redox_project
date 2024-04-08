#for d in /home/hagen/anaconda3/bin/heme_redox_project/hemejpn/database/pdb/from_rcsb_old/*; do
for d in /home/hagen/anaconda3/bin/heme_redox_project/hemejpn/database/pdb/prepared2/*; do
    
    echo "$d"
	cat "$d"| grep "CU " |wc
done

