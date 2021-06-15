#!/bin/bash

# Dim_Bright
samtools merge -n -O BAM -@ 16 ./output/dim_bright.bam ./input/*.bam &&\
    sambamba sort -m 8G --tmpdir /mnt/IM/tmp/ -o ./output/Dim_Bright.bam -t 10 -p output/dim_bright.bam &&\
    rm -rf ./output/dim_bright.bam &&\
    samtools index ./output/Dim_Bright.bam &&\
    samtools view -h ./output/Dim_Bright.bam |\
	awk 'function abs(v) {return v < 0 ? -v : v}; {if(abs($9)>=40 && abs($9)<=140 || $1 ~ /^@/){ print $0; }}' |\
	samtools view -bh -o ./output/Dim_Bright_NFF.bam - &&\
    samtools index ./output/Dim_Bright_NFF.bam

# Dim
samtools merge -n -@ 16 ./output/dim.bam ./input/*Dim*.bam &&\
    sambamba sort -m 8G --tmpdir /mnt/IM/tmp/ -o ./output/Dim.bam -t 10 -p ./output/dim.bam &&\
    rm -rf ./output/dim.bam &&\
    samtools index ./output/Dim.bam &&\
    samtools view -h ./output/Dim.bam |\
	awk 'function abs(v) {return v < 0 ? -v : v}; {if(abs($9)>=40 && abs($9)<=140 || $1 ~ /^@/){ print $0; }}' |\
	samtools view -bh -o ./output/Dim_NFF.bam - &&\
    samtools index ./output/Dim_NFF.bam

# Bright
samtools merge -n -@ 16 ./output/bright.bam ./input/*Bright*.bam &&\
    sambamba sort -m 8G --tmpdir /mnt/IM/tmp/ -o ./output/Bright.bam -t 10 -p ./output/bright.bam &&\
    rm -rf ./output/bright.bam &&\
    samtools index ./output/Bright.bam &&\
    samtools view -h ./output/Bright.bam |\
	awk 'function abs(v) {return v < 0 ? -v : v}; {if(abs($9)>=40 && abs($9)<=140 || $1 ~ /^@/){ print $0; }}' |\
	samtools view -bh -o ./output/Bright_NFF.bam - &&\
    samtools index ./output/Bright_NFF.bam
