#!/bin/bash

## All samples
samtools merge -n -@ 16 --write-index -O BAM ./output/pnd_adult.bam ./output/*.bam &&\
 samtools sort -m 8G -O BAM -@ 16 -o ./output/PND_Adult.bam ./output/pnd_adult.bam &&\
 rm -rf ./output/pnd_adult.ba* &&\
 samtools index ./output/PND_Adult.bam

## PND15 samples
samtools merge -n -@ 16 --write-index -O BAM ./output/PND15.bam ./output/PND*.bam &&\
 samtools sort -m 8G -O BAM -@ 16 -o ./output/PND15.bam ./output/pnd15.bam &&\
 rm -rf ./output/pnd15.ba* &&\
 samtools index ./output/PND15.bam

## Adult samples
samtools merge -n -@ 16 --write-index -O BAM ./output/adult.bam ./output/Adult*.bam &&\
 samtools sort -m 8G -O BAM -@ 16 -o ./output/Adult.bam ./output/adult.bam &&\
 rm -rf ./output/adult.ba* &&\
 samtools index ./output/Adult.bam
