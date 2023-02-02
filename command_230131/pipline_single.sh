#!/bin/bash
source bulk.sh
#用于单端测序
work_dir=$(pwd)
cd ..
all_dir=$(pwd)
cd ${work_dir}
species="mouse"
group_1_name="supwt1" #control
group_2_name="supmu2" #treatment
group1="supwt1"
group2="supmu2"
sample="${group1} ${group2}"
reference_rna_dir="/data1/shenluemo"
reference_dir="/data1/shenluemou/reference/mapping/hisat2/mouse/grcm38/mouse_grcm38"
gtf_dir="/data1/shenluemou/reference/gtf/mouse/grcm38/mouse_grcm38.gtf"
bed_dir="/data1/shenluemou/reference/bed/mouse/grcm38/mouse_grcm38.bed"

mkdir ../result/
#-------fastp---------------
mkdir ../result/0fastp
for name in ${sample};do
	fastp_control ${all_dir}/data ${name} ${all_dir}/result/0fastp
	echo "the fastp of ${name} finished"
done

#-------de_rrna------------
mkdir ../result/1de
for name in ${sample};do
	de_rna ${all_dir}/result/0fastp $name ${all_dir}/result/1de $reference_rna_dir
	echo "the de_rrna of "$name" finished"
done

#-------bw--------------------
mkdir ../result/2map
for name in ${sample};do
	fastqtobw ${all_dir}/result/1de ${name} ${reference_dir} ${all_dir}/result/2map
	echo "the mapping of "${name}" finished"
done

#--------------rnaqc-------------
mkdir ../result/3qc
for name in ${sample};do
	rnaqc ${all_dir}/result/2map ${name} ${all_dir}/result/3qc ${gtf_dir} ${bed_dir}
	echo "the bam-qc of "${name}" finished"
done

#---------------abundance----------
mkdir ../result/4abundance
for name in ${sample};do
	estimate_abundance ${all_dir}/result/2map ${name} ${gtf_dir} ${all_dir}/result/4abundance
	echo "the abundance of "${name}" finished"
done
mkdir ../result/5matrix
total_abundance ${all_dir}/result/4abundance ${all_dir}/result/5matrix
read_matrix ${all_dir}/result/5matrix ${all_dir}/result/4abundance

#----------------saturation------------
saturation ${all_dir}/result/5matrix/qualimeta.txt ${all_dir}/result/5matrix/


