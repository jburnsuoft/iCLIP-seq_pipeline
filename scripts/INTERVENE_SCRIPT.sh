#!/bin/sh

#  INTERVENE_SCRIPT.sh
#  
#
#  Created by james burns on 2021-10-18.
#  


intervene pairwise -i *.top1000.bed --compute frac -o newfolder-frac_top1000/

intervene pairwise -i *.top1000.bed --compute jaccard --sort -o newfolder-jaccard_top1000/
