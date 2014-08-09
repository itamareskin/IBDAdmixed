#!/bin/bash

data_dir=/a/home/cc/cs/itamares/Data/ibdadmxed/New4
pop1_prefix=../New2/HapMap3_CEU_chr2
pop2_prefix=../New2/HapMap3_YRI_chr2
data_prefix=ceu.tsi.yri.lwk.half2
output_name=$data_prefix.1
beagle_dag1=$pop1_prefix.$pop1_prefix.bgl.dag.gz

# run Beagle to create dag model
python ibdadmx.py ped2bgl $data_dir/$pop1_prefix
python ibdadmx.py bglmodel $data_dir/$pop1_prefix

# run GERMLINE
python ibdadmx.py germline $data_dir/$data_prefix.genos

# run ibdadmxed
nohup python ibdadmx.py ibd $data_dir/$data_prefix.genos $data_dir/$output_name $data_dir/$beagle_dag1 -k 1 -a 1 --pairs-file $data_dir/$data_prefix.trueibd.pairs.txt -p 15 --set-ibd-trans 1e-5 1 --germline-file $data_dir/$data_prefix.match -m -50 -w 100 &

# run Naive Model
nohup python ibdadmx.py ibd $data_dir/$data_prefix.genos $data_dir/$output_name $data_dir/$beagle_dag1 -k 1 -a 1 --pairs-file $data_dir/$data_prefix.trueibd.pairs.txt -p 15 --set-ibd-trans 1e-5 1 --germline-file $data_dir/$data_prefix.match -m -50 -w 100 --naive-model &

# run Beagle 3

# run Beagle 4

# run Parente

# run Parente2

# run Germline

# calc stats
python ibdadmx.py stats $data_dir/$data_prefix.trueibd.txt $data_dir/$output_name.ibdadmixed.txt
