# RNASeq_ref_pipe
Run RNASeq analysis with reference genome based on the Tophat-Cufflinks pipeline
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

# 1. input.txt, list the sample name and clean reads in each line, seperated by tab or space
sample1    sample1_R1.fq.gz   sample1_R2.fq.gz
sample2    ample2_R1.fq.gz    sample2_R2.fq.gz
sample3    sample3_R1.fq.gz   sample3_R3.fq.gz

OPTIONS:
    -i input file, list the sample name and clean reads in each line, default input.txt
    -g reference genome, default ref.genome.fasta
    -G gtf/gff file with the reference genome, default ref.genome.gtf
    -t thread, default 10
    -h/? show help of this script
Example:
    Run_RNASeq_ref_pipe_v1.sh -i input.txt -g ref.genome.fasta -G ref.genome.gtf -t 10
