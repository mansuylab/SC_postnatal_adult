#!/bin/bash

#samtools merge -n -@ 8 ./output/pnd_adult.bam ./input/*.bam &&  samtools sort -o ./output/PND_Adult.bam -O BAM -@ 8 ./output/pnd_adult.bam && rm -rf ./output/pnd_adult.bam && samtools index ./output/PND_Adult.bam
samtools merge -n -@ 8 ./output/PND15.bam ./input/PND*.bam &&  samtools sort -o ./output/PND15.bam -O BAM -@ 8 ./output/pnd15.bam && rm -rf ./output/pnd15.bam && samtools index ./output/PND15.bam
samtools merge -n -@ 8 ./output/adult.bam ./input/Adult*.bam &&  samtools sort -o ./output/Adult.bam -O BAM -@ 8 ./output/adult.bam && rm -rf ./output/adult.bam && samtools index ./output/Adult.bam
