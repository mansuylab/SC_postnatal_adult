#!/bin/bash

# PND15_Adult
#samtools merge -n -O BAM -@ 16 ./output/pnd_adult.bam ./input/*.bam &&\
#    sambamba sort -m 8G --tmpdir /mnt/IM/tmp/ -o ./output/PND15_Adult.bam -t 10 -p output/pnd_adult.bam &&\
#    rm -rf ./output/pnd_adult.bam &&\
#    samtools index ./output/PND15_Adult.bam &&\
#    samtools view -h ./output/PND15_Adult.bam |\
#	awk 'function abs(v) {return v < 0 ? -v : v}; {if(abs($9)>=40 && abs($9)<=140 || $1 ~ /^@/){ print $0; }}' |\
#	samtools view -bh -o ./output/PND15_Adult_NF.bam - &&\
#    samtools index ./output/PND15_Adult_NF.bam

# PND15
#samtools merge -n -@ 16 ./output/pnd15.bam ./input/PND15*.bam &&\
#    sambamba sort -m 8G --tmpdir /mnt/IM/tmp/ -o ./output/PND15.bam -t 10 -p ./output/pnd15.bam &&\
#    rm -rf ./output/pnd15.bam &&\
#    samtools index ./output/PND15.bam &&\
#    samtools view -h ./output/PND15.bam |\
#	awk 'function abs(v) {return v < 0 ? -v : v}; {if(abs($9)>=40 && abs($9)<=140 || $1 ~ /^@/){ print $0; }}' |\
#	samtools view -bh -o ./output/PND15_NF.bam - &&\
#    samtools index ./output/PND15_NF.bam

# Adults
#samtools merge -n -@ 16 ./output/adult.bam ./input/Adult*.bam &&\
    sambamba sort -m 8G --tmpdir /mnt/IM/tmp/ -o ./output/Adult.bam -t 10 -p ./output/adult.bam &&\
    rm -rf ./output/adult.bam &&\
    samtools index ./output/Adult.bam &&\
    samtools view -h ./output/Adult.bam |\
	awk 'function abs(v) {return v < 0 ? -v : v}; {if(abs($9)>=40 && abs($9)<=140 || $1 ~ /^@/){ print $0; }}' |\
	samtools view -bh -o ./output/Adult_NF.bam - &&\
    samtools index ./output/Adult_NF.bam
