############################
## Author: Deepak Tanwar  ##
## Date: 20191008	  ##
############################

# This makefile will run the deepTools for filtering bamfile for MAPQ 30.

SHELL:=/bin/bash
source_dir=./input
target_dir=./output

files := $(wildcard $(source_dir)/*.bam)
targets := $(patsubst $(source_dir)/%.bam, $(target_dir)/%.bam, $(files))

all: $(targets)

$(target_dir)/%.bam: $(source_dir)/%.bam
	alignmentSieve --version 2> ./log/$(basename $(notdir $@)).log >> ./log/$(basename $(notdir $@)).log && alignmentSieve --minMappingQuality 30 -p 16 -b $< -o $(basename $@)_tmp.bam 2>> ./log/$(basename $(notdir $@)).log && rm -rf /mnt/IM/tmp/_deep* && samtools sort -m 8G -@ 16 -O BAM -o $@ $(basename $@)_tmp.bam && samtools index $@ && rm -rf $(basename $@)_tmp.bam
