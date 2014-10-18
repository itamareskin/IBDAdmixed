#!/bin/bash

code_dir=/a/home/cc/cs/itamares/Projects/ibdadmixed
data_dir=/home/nasheran/itamares/Data/IBDAdmixed/New5 #/a/home/cc/cs/itamares/Data/IBDAdmixed/New4
pop1_prefix=HapMap3_CEU_chr2
pop2_prefix=HapMap3_YRI_chr2
pop3_prefix=HapMap3_JPT+CHB_chr2
data_prefix=ceu.yri
output_name=$data_prefix.t2
beagle_dag1=$pop1_prefix.$pop1_prefix.bgl.dag.gz
beagle_dag2=$pop2_prefix.$pop2_prefix.bgl.dag.gz
beagle_dag3=$pop3_prefix.$pop3_prefix.bgl.dag.gz
plink=/home/nasheran/itamares/Software/plink-1.07-x86_64/plink
parente=/home/nasheran/itamares/Software/parente

# download HapMap files (run this in data_dir)
#python $code_dir/scripts/loadHapMap3.py 

# convert to plink format
python $code_dir/scripts/saveHapMapPlink.py $data_dir/$pop1_prefix.pop
python $code_dir/scripts/saveHapMapPlink.py $data_dir/$pop2_prefix.pop
python $code_dir/scripts/saveHapMapPlink.py $data_dir/$pop3_prefix.pop

# run Beagle to create dag model
#python $code_dir/scripts/ibdadmx.py ped2bgl $data_dir/$pop1_prefix
python $code_dir/scripts/ibdadmx.py bglmodel $data_dir/$pop1_prefix --scale 4
#python $code_dir/scripts/ibdadmx.py ped2bgl $data_dir/$pop2_prefix
python $code_dir/scripts/ibdadmx.py bglmodel $data_dir/$pop2_prefix --scale 4
#python $code_dir/scripts/ibdadmx.py ped2bgl $data_dir/$pop3_prefix
python $code_dir/scripts/ibdadmx.py bglmodel $data_dir/$pop3_prefix --scale 4

# run simulation
python $code_dir/scripts/simple_simulation.py $data_dir/$pop1_prefix.pop $data_dir/$pop1_prefix.map $data_dir/$data_prefix -a 1 -n 100 -i 80 -e 0.005

# run GERMLINE (necessary for ibdadmx)
python $code_dir/scripts/ibdadmx.py germline $data_dir/$data_prefix.genos $data_dir/$data_prefix.genos --bits 50 --err-hom 2 --err-het 4 --min-m 0.1

# run ibdadmx#python $code_dir/scripts/ibdadmx.py ibd $data_dir/$data_prefix.genos $data_dir/$output_name $data_dir/$beagle_dag1 -k 1 -a 1 --pairs-file $data_dir/$data_prefix.trueibd.pairs.txt -p 1 --set-ibd-trans 1e-5 1 --germline-file $data_dir/$data_prefix.genos.match -m -50 -w 100
nohup python $code_dir/scripts/ibdadmx.py ibd $data_dir/$data_prefix.genos $data_dir/$output_name $data_dir/$beagle_dag1 $data_dir/$beagle_dag2 -k 2 --pairs-file $data_dir/$data_prefix.trueibd.pairs.txt --set-ibd-trans 1e-5 1 1e-5 1 --germline-file $data_dir/$data_prefix.genos.match -m -50 -w 100 -p 16 &

# run Naive Model
nohup python $code_dir/scripts/ibdadmx.py ibd $data_dir/$data_prefix.genos $data_dir/$output_name.naive $data_dir/$beagle_dag1 $data_dir/$beagle_dag2 $data_dir/$beagle_dag3 -k 3 -a 0.2 0.4 0.4 --pairs-file $data_dir/$data_prefix.trueibd.pairs.txt --set-ibd-trans 2e-4 1 1e-5 1 1e-5 1 --germline-file $data_dir/$data_prefix.genos.match -m -50 -w 100 --naive-model --rerun -p 15 &

# run Beagle 3 fastIBD
python $code_dir/scripts/ibdadmx.py ped2bgl $data_dir/$data_prefix.genos
nohup python $code_dir/scripts/ibdadmx.py beagle3 $data_dir/$data_prefix.genos $data_dir/$output_name $data_dir/$pop1_prefix --nruns 7 --fastIBD-threshold 1e-4 &
gunzip $data_dir/$output_name.*$data_prefix.genos.bgl.fibd.gz
awk -v maxscore="1e-3" 'BEGIN{split("",ibd)}$5<=maxscore{split($1,a,"."); split($2,b,"."); h1=int((a[2]-3)/2); h2=int((b[2]-3)/2); p=h1","h2; v=$3-1","$4-1","$5; if (p in ibd) ibd[p]=ibd[p]";"v; else ibd[p]=v}END{for (i in ibd) print i":"ibd[i]}' $data_dir/$output_name.*$data_prefix.genos.bgl.fibd > $data_dir/$output_name.genos.beagle3.ibd.txt

# run Beagle 4
awk 'BEGIN{OFS=" "}NR==1{$1="I"; $2="id"; print; next}{for (i=3;i<=NF;i++) {if ($i=="1") $i="A"; else $i="G"}; print;}' $data_dir/$data_prefix.genos.bgl > $data_dir/$data_prefix.genos.fixed.bgl
awk '{print $2,$4,"A","G"}' $data_dir/$data_prefix.genos.map > $data_dir/$data_prefix.genos.fixed.markers
python $code_dir/scripts/ibdadmx.py bgl2vcf $data_dir/$data_prefix.genos.fixed
sed 's/\//\|/g' $data_dir/$data_prefix.genos.fixed.vcf > $data_dir/$data_prefix.genos.vcf
python $code_dir/scripts/ibdadmx.py beagle4 $data_dir/$data_prefix.genos
awk 'BEGIN{split("",ibd)}NR==FNR{pos[$4]=NR-1; next}{h1=$1; h2=$3; p=h1","h2; v=pos[$6]","pos[$7]","$8; if (p in ibd) ibd[p]=ibd[p]";"v; else ibd[p]=v}END{for (i in ibd) print i":"ibd[i]}' $data_dir/$data_prefix.genos.map $data_dir/$data_prefix.genos.out.ibd > $data_dir/$output_name.genos.beagle4.ibd.txt

# run Parente
$plink --noweb --file $data_dir/$pop1_prefix --recode --transpose --out $data_dir/$pop1_prefix
$plink --tfile $data_dir/$pop1_prefix --noweb --freq --out $data_dir/$pop1_prefix
$parente tool tped_to_ints $data_dir/$pop1_prefix $data_dir/$pop1_prefix.frq hap $data_dir/$pop1_prefix
$plink --noweb --file  $data_dir/$data_prefix.genos --recode --transpose --out $data_dir/$data_prefix
$plink --tfile $data_dir/$data_prefix --noweb --freq --out $data_dir/$data_prefix
$parente tool tped_to_ints $data_dir/$data_prefix $data_dir/$data_prefix.frq geno $data_dir/$data_prefix
$parente train -w 20 -t 16 $data_dir/$pop1_prefix.hap $data_dir/$data_prefix.tped $data_dir/$data_prefix
$parente infer -s 2 -t 16 $data_dir/$data_prefix $data_dir/$data_prefix.geno -20 $data_dir/$data_prefix.2
awk 'NR==FNR{a[$3]=FNR; next;}{print $1,a[$2],a[$3]}' $data_dir/$pop1_prefix.tped $data_dir/$data_prefix.2.fblock | awk 'NR==FNR{a[$1]=$2;b[$1]=$3;next}{print $2,$3,a[$4],b[$4],$5}' - $data_dir/$data_prefix.2.ibd | awk 'BEGIN{split("",ibd)}{h1=$1; h2=$2; p=h1","h2; v=$3","$4","$5; if (p in ibd) ibd[p]=ibd[p]";"v; else ibd[p]=v}END{for (i in ibd) print i":"ibd[i]}' > $data_dir/$data_prefix.parente.ibd.txt

# run Parente2

# run Germline

for i in {21..50}
do
	min_m=$(awk -v i=$i 'BEGIN { print ((i) / 10) }')
	#echo $min_m
	python $code_dir/scripts/ibdadmx.py germline $data_dir/$data_prefix.genos $data_dir/$output_name.$i --bits 64 --min-m $min_m --err-hom 0 --err-het 0
done

## calc stats
# ibdadmx
python $code_dir/scripts/ibdadmx.py stats $data_dir/$pop1_prefix.map $data_dir/$data_prefix.trueibd.txt $data_dir/$output_name.ibdadmixed.txt --lod-score --min-score-point -50 --max-score-point 200 --num-score-points 150
python $code_dir/scripts/ibdadmx.py stats $data_dir/$pop1_prefix.map $data_dir/$data_prefix.trueibd.txt $data_dir/$output_name.ibdadmixed.0.5.txt --lod-score --min-score-point -50 --max-score-point 200 --filter-est-length 0.5 --num-score-points 100
# naive
python $code_dir/scripts/ibdadmx.py stats $data_dir/$pop1_prefix.map $data_dir/$data_prefix.trueibd.txt $data_dir/$output_name.naive.ibdadmixed.txt --lod-score --min-score-point -50 --max-score-point 200
# beagle3
python $code_dir/scripts/ibdadmx.py stats $data_dir/$pop1_prefix.map $data_dir/$data_prefix.trueibd.txt $data_dir/$output_name.genos.beagle3.ibd.txt --beagle --min-score-point " -50" --max-score-point " -3" --compare-same-inds --min-score -10000 --max-score 0
# beagle4
python $code_dir/scripts/ibdadmx.py stats $data_dir/$pop1_prefix.map $data_dir/$data_prefix.trueibd.txt $data_dir/$output_name.genos.beagle4.ibd.txt --lod-score
# parente
python $code_dir/scripts/ibdadmx.py stats $data_dir/$pop1_prefix.map $data_dir/$data_prefix.trueibd.txt $data_dir/$data_prefix.parente.ibd.txt --parente --lod-score --min-score -20 --max-score 50
# germline
python $code_dir/scripts/ibdadmx.py stats $data_dir/$pop1_prefix.map $data_dir/$data_prefix.trueibd.txt $data_dir/$output_name.germline.ibd.txt
for i in {21..50}
do
	min_m=$(awk -v i=$i 'BEGIN { print ((i) / 10) }')
	#echo $min_m
	python $code_dir/scripts/ibdadmx.py stats $data_dir/$pop1_prefix.map $data_dir/$data_prefix.trueibd.txt $data_dir/$output_name.$i.germline.ibd.txt --one-line
done

awk '{print NR,$1,$2,$3}' `ls -v $data_dir/$output_name.*.germline.ibd.txt.stats.txt` > $data_dir/$output_name.germline.ibd.txt.stats.txt



for l in 2
do
	let r=$l-1
	# ibdadmx
	python $code_dir/scripts/ibdadmx.py stats $data_dir/$pop1_prefix.map $data_dir/$data_prefix.trueibd.txt $data_dir/$output_name.ibdadmixed.txt --lod-score --min-score-point " -50" --max-score-point 200 --num-score-points 75 --min-length $r.5 --max-length $l.5 --out $data_dir/$output_name.ibdadmixed.txt.$l.cM
	# naive
	python $code_dir/scripts/ibdadmx.py stats $data_dir/$pop1_prefix.map $data_dir/$data_prefix.trueibd.txt $data_dir/$output_name.naive.ibdadmixed.txt --lod-score --min-score-point " -50" --max-score-point 200 --num-score-points 75 --min-length $r.5 --max-length $l.5 --out $data_dir/$output_name.naive.ibdadmixed.txt.$l.cM
	# beagle3
	python $code_dir/scripts/ibdadmx.py stats $data_dir/$pop1_prefix.map $data_dir/$data_prefix.trueibd.txt $data_dir/$output_name.genos.beagle3.ibd.txt --beagle --min-score-point " -50" --max-score-point " -3" --compare-same-inds --min-score -10000 --max-score 0 --min-length $r.5 --max-length $l.5 --out $data_dir/$output_name.genos.beagle3.ibd.txt.$l.cM
	# parente
	python $code_dir/scripts/ibdadmx.py stats $data_dir/$pop1_prefix.map $data_dir/$data_prefix.trueibd.txt $data_dir/$data_prefix.parente.ibd.txt --parente --lod-score --min-score-point " -20" --max-score-point 50 --num-score-points 75 --min-score -20 --max-score 50 --min-length $r.5 --max-length $l.5 --out $data_dir/$data_prefix.parente.ibd.txt.$l.cM
# # germline
# for i in {1..50}
# do
# 	min_m=$(awk -v i=$i 'BEGIN { print ((i) / 10) }')
# 	#echo $min_m
# 	python $code_dir/scripts/ibdadmx.py stats $data_dir/$pop1_prefix.map $data_dir/$data_prefix.trueibd.txt $data_dir/$output_name.$i.germline.ibd.txt --one-line --min-length 1.5 --max-length 2.5 --out $data_dir/$output_name.$i.germline.ibd.txt.2cM
# done
# awk '{print NR,$1,$2,$3}' `ls -v $data_dir/$output_name.*.germline.ibd.txt.stats.txt.2cM` > $data_dir/$output_name.germline.ibd.txt.2cM.stats.txt
done 