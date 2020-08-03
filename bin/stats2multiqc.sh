#!/bin/bash


function usage {                                                                                                                                                                                                echo -e "usage : stats2multiqc.sh -s SAMPLE_PLAN -a ALIGNER -d STRANDNESS -f PATERNAL -m MATERNAL [-p][-n][-h]" 
    echo -e "Use option -h|--help for more information"
}

function help {
    usage;
    echo
    echo "stat2multiqc.sh"
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -s SAMPLE_PLAN"
    echo "   -a ALIGNER"
    echo "   -d STRANDNESS"
    echo "   -f PATERNAL GENOTYPE"
    echo "   -m MATERNAL GENOTYPE"
    echo "   [-p] paired-end"
    echo "   [-n] N-mask"
    echo "   [-h]: help"
    exit;
}

nmask=0
is_pe=0
while getopts "s:a:d:f:m:nph" OPT
do
    case $OPT in
        s) splan=$OPTARG;;
        a) aligner=$OPTARG;;
        d) strandness=$OPTARG;;
	f) paternal=$OPTARG;;
	m) maternal=$OPTARG;;
	p) is_pe=1;;
	n) nmask=1;;
        h) help ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            usage
            exit 1
            ;;
    esac
done

if  [[ -z $splan ]]; then
    usage
    exit
fi

echo -e "Sample_id,Sample_name,Genotypes,Number_of_reads,Strandness,Number_of_aligned_reads,Percent_of_aligned_reads,Number_paternal,Percent_paternal,Number_maternal,Percent_maternal,Number_unassigned,Percent_unassigned,Number_conflicting,Percent_conflicting" > mq.stats

if [ ${nmask} == "1" ]; then
    genomeBase=${paternal}_${maternal}_nmask_genome
else
    genomeBase=${paternal}_${maternal}
fi

## Catch sample names
all_samples=$(awk -F, '{print $1}' $splan)

echo -e "Sample_id,Sample_name,Genotypes,Number_of_reads,Strandness,Number_of_aligned_reads,Percent_of_aligned_reads,Number_paternal,Percent_paternal,Number_maternal,Percent_maternal,Number_unassigned,Percent_unassigned,Number_conflicting,Percent_conflicting" > mq.stats

for sample in $all_samples
do
    ##id
    sname=$(awk -F, -v sname=$sample '$1==sname{print $2}' $splan)

    ## N-mask mapping
    if [ ${nmask} == "1" ]; then
    
	if [ $aligner == "star" ]; then
	    n_reads=$(grep "Number of input reads" alignment/${sample}_${genomeBase}Log.final.out | cut -d"|" -f 2 | sed -e 's/\t//g')
	elif [ $aligner == "tophat2" ]; then
	    n_reads=$(grep "Input" alignment/${sample}_${genomeBase}.align_summary.txt | uniq | cut -d: -f2 | sed -e 's/ //g')
	elif [ $aligner == "bowtie2" ]; then
	    n_reads=$(grep "reads;" alignment/${sample}_${genomeBase}_bowtie2.log | awk '{print $1}')
	elif [[ $aligner == "hisat2" && $is_pe == "1" ]]; then
	    n_reads=$(grep "Total pairs" alignment/${sample}.hisat2_summary.txt | cut -d: -f2 | sed -e 's/ //g')
	elif [[ $aligner == "hisat2" && $is_pe == "0" ]]; then
	    n_reads=$(grep "Total reads" alignment/${sample}.hisat2_summary.txt | cut -d: -f2 | sed -e 's/ //g')
	fi

	if [ $aligner == "tophat2" ]; then
	    n_mapped=$(grep "Aligned pairs" alignment/${sample}_${genomeBase}.align_summary.txt | cut -d: -f 2 | sed -e 's/ //g')
	    n_multi=$(grep -a2 "Aligned pairs" alignment/${sample}_${genomeBase}.align_summary.txt | grep "multiple" | awk -F" " '{print $3}')
	    n_unique=$(($n_mapped - $n_multi))
	elif [ $aligner == "star" ]; then
	    n_unique=$(grep "Uniquely mapped reads number" alignment/${sample}_${genomeBase}Log.final.out | cut -d"|" -f 2 | sed -e 's/\t//g')
	    n_multi=$(grep "Number of reads mapped to multiple loci" alignment/${sample}_${genomeBase}Log.final.out | cut -d"|" -f 2 | sed -e 's/\t//g')
	    n_mapped=$(($n_unique + $n_multi))
	elif [ $aligner == "hisat2" ]; then
	    n_unique=$(grep " 1 time" alignment/${sample}_${genomeBase}.hisat2_summary.txt | cut -d: -f 2 | sed -e 's/ //g' | awk -F"(" 'BEGIN{s=0}{s=s+$1}END{print s}')
	    n_multi=$(grep ">1 time" alignment/${sample}_${genomeBase}.hisat2_summary.txt | cut -d: -f 2 | sed -e 's/ //g' | awk -F"(" 'BEGIN{s=0}{s=s+$1}END{print s}')
	    n_mapped=$(($n_unique + $n_multi))
	elif [ $aligner == "bowtie2" ]; then
	    if [[ $is_pe == "1" ]]; then
		n_uniq=$(grep "concordantly exactly" alignment/${sample}_${genomeBase}_bowtie2.log | awk '{print $1}')
		n_multi=$(grep "concordantly >1" alignment/${sample}_${genomeBase}_bowtie2.log | awk '{print $1}')
	    else
		n_uniq=$(grep "exactly" alignment/${sample}_${genomeBase}_bowtie2.log | awk '{print $1}')
		n_multi=$(grep ">1" alignment/${sample}_${genomeBase}_bowtie2.log | awk '{print $1}')
	    fi
	    n_mapped=$(($n_uniq + $n_multi))
	else
	    n_mapped='NA'
	fi
	
    ## Parental Mapping
    else
	n_reads=$(grep "Total" tag/${sample}_mergeAlignReport.log | awk -F"\t" '{print $2}' | awk -F\( '{print $1}')
	n_mapped=$(grep "Reads mapped:" tag/${sample}_mergeAlignReport.log | awk -F"\t" '{print $2}' | awk -F\( '{print $1}')
	if [[ ${is_pe} == "1" ]]; then
	    n_reads=$((n_reads/2))
	    n_mapped=$((n_mapped/2))
	fi
    fi

    ## AS
    n_un=$(grep "Unassigned" snpsplit/${sample}*sort.mqc | awk -F, '{print $2}' | sed -e 's/ //g')
    n_g1=$(grep "Genome 1" snpsplit/${sample}*sort.mqc | awk -F, '{print $2}' | sed -e 's/ //g')
    n_g2=$(grep "Genome 2" snpsplit/${sample}*sort.mqc | awk -F, '{print $2}' | sed -e 's/ //g')
    n_conf=$(grep "Conflicting" snpsplit/${sample}*sort.mqc | awk -F, '{print $2}' | sed -e 's/ //g')

    ## Calculate percentage
    p_mapped=$(echo "${n_mapped} ${n_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    p_un=$(echo "${n_un} ${n_mapped}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    p_g1=$(echo "${n_g1} ${n_mapped}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    p_g2=$(echo "${n_g2} ${n_mapped}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    p_conf=$(echo "${n_conf} ${n_mapped}" | awk ' { printf "%.*f",2,$1*100/$2 } ')

    echo -e ${sample},${sname},${paternal}"/"${maternal},${n_reads},${strandness},${n_mapped},${p_mapped},${n_g1},${p_g1},${n_g2},${p_g2},${n_un},${p_un},${n_conf},${p_conf} >> mq.stats
done


## Heatmap allelic ratio
if [ -d asratio ]; then
    header=0
    for i in asratio/* 
    do 
	if [ ${header} == 0 ]
	then
	    awk -F"," '{printf ","$1}END{printf "\n"}' $i > asratio_heatmap.mqc
	    header=1
	fi
	id=$(basename $i | sed -e 's/_average.mqc//')
	awk -F, -v id=${id} 'BEGIN{printf id}{printf ","$2}END{printf "\n"}' $i >> asratio_heatmap.mqc
    done

fi
