###########################
## Author: Deepak Tanwar ##
## Date: 20210111        ##
###########################

## This makefile will mark duplicated reads of processed bam file.

source_dir=./input
target_dir=./output

files := $(wildcard $(source_dir)/*.bam)
targets := $(patsubst $(source_dir)/%.bam, $(target_dir)/%.bam, $(files))

all: $(targets)

$(target_dir)/%.bam: $(source_dir)/%.bam
	sambamba markdup -t 8 -p --tmpdir /mnt/IM/tmp/ $< $@ 2>./log/$(basename $(notdir $@)).log && samtools index $@
