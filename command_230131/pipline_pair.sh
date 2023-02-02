#!/bin/bash
#用于双端测序
source bulk.sh
source bulk_or_pair.sh

work_dir=$(pwd)
cd ..
all_dir=$(pwd)
cd ${work_dir}
species="mouse"
group_1_name="g1"  #control
group_2_name="g2"	  #treatment
group1="zygote_g1_3"
group2="zygote_g2_3"
sample="${group1} ${group2}"
reference_rna_dir="/data1/shenluemou/reference/mapping/hisat2/mouse/grcm38/mouse_grcm38_rRNA"
reference_dir="/data1/shenluemou/reference/mapping/hisat2/mouse/grcm38/mouse_grcm38"
gtf_dir="/data1/shenluemou/reference/gtf/mouse/grcm38/mouse_grcm38.gtf"
bed_dir="/data1/shenluemou/reference/bed/mouse/grcm38/mouse_grcm38.bed"
nsample="multi" #multi or single
mkdir ../result/

#-------fastp---------------
mkdir ../result/fastp_1
for name in ${sample};do
	bop_fastp_control ${all_dir}/data ${name} ${all_dir}/result/fastp_1
	echo "the fastp of "${name}" finished"
done

#------de_rrna-------------
mkdir ../result/de_2
for name in ${sample};do
	bop_derna ${all_dir}/result/fastp_1 $name $reference_rna_dir ${all_dir}/result/de_2
	echo "the derrna of "$name" finished"
done

#-------bw--------------------
mkdir ../result/map_2
for name in ${sample};do
	bop_map ${all_dir}/result/de_2 ${name} ${reference_dir} ${all_dir}/result/map_2
	echo "the mapping of "${name}" finished"
done

#--------------rnaqc-------------
mkdir ../result/bamqc_3
for name in ${sample};do
	rnaqc ${all_dir}/result/map_2 ${name} ${all_dir}/result/bamqc_3 ${gtf_dir} ${bed_dir}
	echo "the bam-qc of "${name}" finished"
done

#---------------abundance----------
mkdir ../result/abundance_4
for name in ${sample};do
	estimate_abundance ${all_dir}/result/map_2 ${name} ${gtf_dir} ${all_dir}/result/abundance_4
	echo "the abundance of "${name}" finished"
done
mkdir ../result/matrix_5
total_abundance ${all_dir}/result/abundance_4 ${all_dir}/result/matrix_5
read_matrix ${all_dir}/result/matrix_5 ${all_dir}/result/abundance_4

#----------------saturation------------
saturation ${all_dir}/result/matrix_5/qualimeta.txt ${all_dir}/result/matrix_5/

#----------------diff_ana------------------
mkdir ../result/diff_6
touch ${all_dir}/result/matrix_5/matrix.csv
echo "Sample,Group" > ${all_dir}/result/matrix_5/matrix.csv
for name in ${group1};do
        echo "${group_1_name},${name}" >> ${all_dir}/result/matrix_5/matrix.csv
done
for name in ${group2};do
        echo "${group_2_name},${name}" >> ${all_dir}/result/matrix_5/matrix.csv
done
if [ $nsample = "multi" ];then
	Rscript different_multi.R -m ${all_dir}/result/matrix_5/matrix.csv -f ${all_dir}/result/matrix_5/gene_fpkm_matrix.csv -c ${all_dir}/result/matrix_5/gene_count_matrix.csv -o ${all_dir}/result/diff_6/ -s $species
fi
if [ $nsample = "single" ];then
	Rscript different_single.R -m ${all_dir}/result/matrix_5/matrix.csv -f ${all_dir}/result/matrix_5/gene_fpkm_matrix.csv -c ${all_dir}/result/matrix_5/gene_count_matrix.csv -o ${all_dir}/result/diff_6/ -s $species
fi


