#!/bin/bash
set -e
set -u

# 设置程序参数的缺省值
# Default parameter
INPUT=input.txt
GENOME=ref.genome.fasta
GTF=ref.genome.gtf
THREAD=10

# 程序功能描述，每次必改程序名、版本、时间等；版本更新要记录清楚，一般可用-h/-?来显示这部分
# Function for script description and usage
usage()
{
cat <<EOF >&2
Usage:
---------------------------------------------------------------------------------------------------
Filename:    Run_RNASeq_ref_pipe_v1.sh
Revision:    1.0
Date:        2018/10/24
Author:      Wei Dong
Email:       1369852697@qq.com
GitHub:      https://github.com/Davey1220/
Description: This script was designed to conduct the RNASeq analysis with reference genome based on the Tophat-Cufflinks pipeline
---------------------------------------------------------------------------------------------------
Copyright:   2018 (c) Wei Dong
License:     GPL
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License \
as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version. \
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty \
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. \
If any changes are made to this script, please mail me a copy of the changes
---------------------------------------------------------------------------------------------------
Version 1.0 2018/10/24

# 输入文件格式和示例
# 1. input.txt, list the sample name and clean reads in each line, seperated by tab or space
sample1    sample1_R1.fq.gz   sample1_R2.fq.gz
sample2    sample2_R1.fq.gz   sample2_R2.fq.gz
sample3    sample3_R1.fq.gz   sample3_R3.fq.gz

OPTIONS:
    -i input file, list the sample name and clean reads in each line, defult input.txt
    -g reference genome, default ref.genome.fasta
    -G gtf/gff file with the reference genome, default ref.genome.gtf
    -t thread, default 10
    -h/? show help of this script
Example:
    Run_RNASeq_ref_pipe_v1.sh -i input.txt -g ref.genome.fasta -G ref.genome.gtf -t 10
EOF
}

# 解释命令行参数，其实是调用了perl语言的getopts包，
# Analysis parameter
while getopts "g:h:i:G:t:" OPTION
do
    case $OPTION in
        g)
            GENOME=$OPTARG
            ;;
        h)
            usage
            exit 1
            ;;
        i)
            INPUT=$OPTARG
            ;;
        G)
            GTF=$OPTARG
            ;;
        t)
            THREAD=$OPTARG
            ;;                        
        ?)
            usage
            exit 1
            ;;
    esac

done

if [ $# -eq 0 ];then
    usage
    exit 1
fi

# check the dependent software
check_software(){
    software=$1
    path=`which $software`
    tmp=`echo $?`
    if [ $tmp == 0 ];then
        echo -e "\033[32m The $software has been installed and it's path is: $path \033[0m"
    else
        echo -e "\033[32m Can't locate the $software, please check and install it first! \033[0m"
    fi
}

#############################################################
#  Step_1: Lets's first make index for the reference genome
#############################################################
bowtie2_index(){
    INDEX_DIR="01bowtie2_index"
    if [ ! -d $INDEX_DIR ];then
        mkdir $INDEX_DIR
    else
        echo "The $INDEX_DIR exist!"
    fi
    check_software bowtie2-build
    # index
    echo `date` "Start indexing the reference genome"
    bowtie2-build $GENOME $INDEX_DIR/${GENOME%%.fasta} >bowtie2_index.log
    sleep 1
    ln -s `pwd`/$GENOME $INDEX_DIR/${GENOME%%sta}
    echo `date` "Indexing done..."
}

#############################################################
#  Step_2: Mapping the clean data on the reference genome
#############################################################
tophat_mapping(){
    Mapping_DIR="02tophat_mapping"
    if [ ! -d $Mapping_DIR ];then
        mkdir $Mapping_DIR
    else
        echo "The $Mapping_DIR exist!"
    fi
    check_software tophat2
    # tophat
    echo `date` "Start mapping the clean data on the reference genome"
    SampleNum=`wc -l $INPUT`
    cat $INPUT|while read line;do
        SampleName=`echo $line|awk '{print $1}'`
        read1=`echo $line|awk '{print $2}'`
        read2=`echo $line|awk '{print $3}'`
        echo $SampleName
        echo $read1
        echo $read2
        tophat2 -G $GTF -p $THREAD -o $Mapping_DIR/tophat_$SampleName $INDEX_DIR/${GENOME%%.fasta} $read1 $read2 2>tophat_${SampleName}.log
        sleep 1
    done
    echo `date` "Tophat mapping done..."
    # check results
    i=`ls -l $Mapping_DIR/tophat_*/accepted_hits.bam|wc -l`
    if [ $i -ne $SampleNum ];then 
        exit 
    else 
        echo "Tophat was done good" 
    fi
}

#############################################################
#  Step_3: Reconstruct the transcripts
#############################################################
cufflinks_assembly(){
    Assembly_DIR="03cufflinks_assembly"
    if [ ! -d $Assembly_DIR ];then
        mkdir $Assembly_DIR
    else
        echo "The $Assembly_DIR exist!"
    fi
    check_software cufflinks
    # cufflinks
    echo `date` "Start reconstruct each transcripts" 
    SampleNum=`wc -l $INPUT`
    cat $INPUT|while read line;do
        SampleName=`echo $line|awk '{print $1}'`
        cufflinks -p $THREAD -g $GTF -o $Assembly_DIR/cufflinks_$SampleName 02tophat_mapping/tophat_$SampleName/accepted_hits.bam 2>cufflinks_${SampleName}.log
        sleep 1
    done
    echo `date` "Cufflinks assembly done..."
    # check results
    i=`ls -l $Assembly_DIR/cufflinks_*/transcripts.gtf|wc -l`
    if [ $i -ne $SampleNum ];then 
        exit 
    else 
        echo "Cufflinks was done good"
    fi
    # cuffmerge
    rm assembly_GTF_list.txt
    cat $INPUT|while read line;do
        SampleName=`echo $line|awk '{print $1}'`
        echo "$Assembly_DIR/cufflinks_$SampleName/transcripts.gtf" >> assembly_GTF_list.txt
    done 
    cuffmerge -g $GTF -s $GENOME -p $THREAD assembly_GTF_list.txt
    echo `date` "Cuffmerge done..."
}

#############################################################
#  Step_4: Quantify the transcripts
#############################################################
cuffquant_quantify(){
    Quantify_DIR="04quantify"
    if [ ! -d $Quantify_DIR ];then
        mkdir $Quantify_DIR
    else
        echo "The $Quantify_DIR exist!"
    fi
    check_software cuffquant
    # cufflinks
    echo `date` "Start quantify each transcripts" 
    SampleNum=`wc -l $INPUT`
    cat $INPUT|while read line;do
        SampleName=`echo $line|awk '{print $1}'`
        cuffquant merged_asm/merged.gtf -p $THREAD 02tophat_mapping/tophat_$SampleName/accepted_hits.bam -o $Quantify_DIR/cuffquant_$SampleName 2>cuffquant_${SampleName}.log
        sleep 1
    done
    echo `date` "Cuffquant quantify done..."
    # check results
    i=`ls -l $Quantify_DIR/cuffquant_*/abundances.cxb|wc -l`
    if [ $i -ne $SampleNum ];then 
        exit
    else 
        echo "Cuffquant was done good"
    fi
}

#############################################################
#  Step_5: Get read count of every gene
#############################################################
get_readcount(){
    Readcount_DIR="05readcount"
    if [ ! -d $Readcount_DIR ];then
        mkdir $Readcount_DIR
    else
        echo "The $Readcount_DIR exist!"
    fi
    check_software htseq-count
    # htseq-count
    echo `date` "Start get read count of every gene"
    SampleNum=`wc -l $INPUT`
    cat $INPUT|while read line;do
        SampleName=`echo $line|awk '{print $1}'`
        samtools view 02tophat_mapping/tophat_$SampleName/accepted_hits.bam >02tophat_mapping/tophat_$SampleName/accepted_hits.sam
        htseq-count 02tophat_mapping/tophat_$SampleName/accepted_hits.sam $GTF >$Readcount_DIR/${SampleName}.gene.count 2>htseq_${SampleName}.count.log
        rm 02tophat_mapping/tophat_$SampleName/accepted_hits.sam
        sleep 1
    done
    echo `date` "htseq-count done..."
    # check results
    i=`ls -l $Readcount_DIR/*gene.count|wc -l`
    if [ $i -ne $SampleNum ];then 
        exit 
    else 
        echo "HTSeq-count was done good"
    fi 
}

##### run the pipeline #####
echo "#############################################################################################"
echo -e "\033[32m##### Step_1: Lets's first make index for the reference genome #####\033[0m"
echo "#############################################################################################"
bowtie2_index

echo "#############################################################################################"
echo -e "\033[32m##### Step_2: Mapping the clean data on the reference genome #####\033[0m"
echo "#############################################################################################"
tophat_mapping

echo "#############################################################################################"
echo -e "\033[32m##### Step_3: Reconstruct the transcripts #####\033[0m]"
echo "#############################################################################################"
cufflinks_assembly

echo "#############################################################################################"
echo -e "\033[32m##### Step_4: Quantify the transcripts #####\033[0m]"
echo "#############################################################################################"
cuffquant_quantify

echo "#############################################################################################"
echo -e "\033[32m##### Step_5: Get read count of every gene #####\033[0m]"
echo "#############################################################################################"
get_readcount

echo "#############################################################################################"
echo "##### Congratulations, all tasks were finished, Run RNASeq_ref_pipe was done... #####"
echo "#############################################################################################"
