#!/bin/bash

# PND_Adult
#samtools merge -n -O BAM -@ 8 ./output/pnd_adult.bam ./input/*.bam &&  samtools sort -o ./output/PND_Adult.bam -O BAM -@ 8 ./output/pnd_adult.bam && rm -rf ./output/pnd_adult.bam && samtools index ./output/PND_Adult.bam
#samtools view -h ./output/PND_Adult.bam | awk 'function abs(v) {return v < 0 ? -v : v}; {if(abs($9)>=40 && abs($9)<=140 || $1 ~ /^@/){ print $0; }}' | samtools view -bh -o ./output/PND_Adult_NF.bam - && samtools index ./output/PND_Adult_NF.bam

# PND15
samtools merge -n -@ 16 ./output/pnd15.bam ./input/PND15*.bam && samtools sort -@ 16 ./output/pnd15.bam ./output/PND15.bam && rm -rf ./output/pnd15.bam && samtools index ./output/PND15.bam
 samtools view -h ./output/PND15.bam | awk 'function abs(v) {return v < 0 ? -v : v}; {if(abs($9)>=40 && abs($9)<=140 || $1 ~ /^@/){ print $0; }}' > ./output/pnd15.bam &&  samtools view -bhS -o ./output/PND15_NF.bam ./output/pnd15.bam && samtools index ./output/PND15_NF.bam && rm -rf ./output/pnd15.bam

# Adults
samtools merge -n -@ 16 ./output/adult.bam ./input/Adult*.bam && samtools sort -@ 16 ./output/adult15.bam ./output/Adult.bam && rm -rf ./output/adult.bam && samtools index ./output/Adult.bam
samtools view -h ./output/Adult.bam | awk 'function abs(v) {return v < 0 ? -v : v}; {if(abs($9)>=40 && abs($9)<=140 || $1 ~ /^@/){ print $0; }}' > ./output/adult.bam &&  samtools view -bhS -o ./output/Adult_NF.bam ./output/adult.bam && samtools index ./output/Adult_NF.bam && rm -rf ./output/adult.bam
