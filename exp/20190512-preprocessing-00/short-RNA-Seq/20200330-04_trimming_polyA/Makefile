############################
## Author: Deepak Tanwar  ##
## Date: 20200330	  ##
############################

# This makefile will run the TrimGalore software to improve the quality of FastQ files

SHELL:=/bin/bash
source_dir=./input
target_dir=./output

files := $(wildcard $(source_dir)/*.fq.gz)
targets := $(patsubst $(source_dir)/%.fq.gz, $(target_dir)/%_trimmedPolyA.fq.gz, $(files))

all: $(targets)

$(target_dir)/%_trimmedPolyA.fq.gz: $(source_dir)/%.fq.gz
#	trim_galore $< -q 30 --length 15 -j 16 -o $(target_dir) > ./log/$(basename $(basename $(notdir $@))).log 2>>./log/$(basename $(basename $(notdir $@))).log
	cutadapt -m 16 -q 30 -u 3 -a AAAAAAAAAA -j 8 --trim-n -o $@ $< > ./log/$(basename $(basename $(notdir $@))).log 2>>./log/$(basename $(basename $(notdir $@))).log
