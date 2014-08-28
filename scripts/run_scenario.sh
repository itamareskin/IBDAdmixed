#!/bin/bash

code_dir=~/Projects/ibdadmixed
data_dir=~/Data/IBDAdmixed #/a/home/cc/cs/itamares/Data/IBDAdmixed/New4
pop1_prefix=HapMap3_CEU_chr2
pop2_prefix=HapMap3_YRI_chr2
data_prefix=ceu.yri
output_name=$data_prefix.1
beagle_dag1=$pop1_prefix.$pop1_prefix.bgl.dag.gz

# download HapMap files (run this in data_dir)
# python $code_dir/scripts/loadHapMap3.py 

# convert to plink format
python $code_dir/scripts/saveHapMapPlink.py $data_dir/$pop1_prefix.pop
python $code_dir/scripts/saveHapMapPlink.py $data_dir/$pop2_prefix.pop

# run Beagle to create dag model
python $code_dir/scripts/ibdadmx.py ped2bgl $data_dir/$pop1_prefix
python $code_dir/scripts/ibdadmx.py bglmodel $data_dir/$pop1_prefix
python $code_dir/scripts/ibdadmx.py ped2bgl $data_dir/$pop2_prefix
python $code_dir/scripts/ibdadmx.py bglmodel $data_dir/$pop2_prefix

# run simulation
python $code_dir/scripts/simple_simulation.py $pop1_prefix.pop $pop2_prefix.pop $pop1_prefix.map $data_prefix -a 0.2 0.8 -n 100 -i 80 -e 0.005

# run GERMLINE
python $code_dir/scripts/ibdadmx.py germline $data_dir/$data_prefix.genos

# run ibdadmxed
python $code_dir/scripts/ibdadmx.py ibd $data_dir/$data_prefix.genos $data_dir/$output_name $data_dir/$beagle_dag1 -k 1 -a 1 --pairs-file $data_dir/$data_prefix.trueibd.pairs.txt -p 1 --set-ibd-trans 1e-5 1 --germline-file $data_dir/$data_prefix.genos.match -m -50 -w 100

# run Naive Model
nohup python $code_dir/scripts/ibdadmx.py ibd $data_dir/$data_prefix.genos $data_dir/$output_name $data_dir/$beagle_dag1 -k 1 -a 1 --pairs-file $data_dir/$data_prefix.trueibd.pairs.txt -p 15 --set-ibd-trans 1e-5 1 --germline-file $data_dir/$data_prefix.genos.match -m -50 -w 100 --naive-model &

# run Beagle 3

# run Beagle 4

# run Parente

# run Parente2

# run Germline

# calc stats
python ibdadmx.py stats $data_dir/$data_prefix.trueibd.txt $data_dir/$output_name.ibdadmixed.txt
