###########################
## Author: Deepak Tanwar ##
## Date:   20190816      ##
###########################

## This makefile will perform quantification of data using Salmon.

SHELL:=/bin/bash
source_dir=./input
target_dir=./output

files := $(wildcard $(source_dir)/*_trimmed.fq.gz)
targets := $(patsubst $(source_dir)/%_trimmed.fq.gz, $(target_dir)/%.bam, $(files))

all: $(targets)

$(target_dir)/%.bam: $(source_dir)/%_trimmed.fq.gz

	salmon -v 2> ./log/$(basename $(notdir $@)).log >> ./log/$(basename $(notdir $@)).log && salmon quant -i ./input/gencode.vM18_plus_repeatsTx_plus_piRNAprecursors/ -l A -r <(zcat $<) -p 8 --seqBias --gcBias --dumpEq -o $(basename $@) 2>> ./log/$(basename $(notdir $@)).log >> ./log/$(basename $(notdir $@)).log
