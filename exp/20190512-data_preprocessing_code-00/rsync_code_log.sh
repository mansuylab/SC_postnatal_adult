#!/bin/bash

rsync -avzC\
      --include="*/"\
      --include='*.R'\
      --include='*.Rmd'\
      --include='Makefile'\
      --include='*.sh'\
      --include='*.html'\
      --include='log'\
      --include='*.log'\
      --exclude="*"\
      /mnt/IM/data/SC_controls/data/20190512-preprocessing-00/\
      ./
