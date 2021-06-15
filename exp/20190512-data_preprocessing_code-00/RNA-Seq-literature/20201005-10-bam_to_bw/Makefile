############################
## Author: Deepak Tanwar  ##
## Date: 20201005	  ##
############################

# Convert bam to bw

SHELL:=/bin/bash
source_dir=./input
target_dir=./output

files := $(wildcard $(source_dir)/*.bam)
targets := $(patsubst $(source_dir)/%.bam, $(target_dir)/%.bw, $(files))

all: $(targets)

$(target_dir)/%.bw: $(source_dir)/%.bam
	bamCoverage -b $< -o $@ -p 16 --normalizeUsing RPKM > ./log/$(basename $(notdir $@)).log 2>>./log/$(basename $(notdir $@)).log
