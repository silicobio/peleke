#!/bin/bash

## Script to run Humatch on the generated antibody sequences

# git clone https://github.com/oxpig/Humatch.git
# cd ./Humatch && pip install .


Humatch-classify -i generated_antigen_seqs.csv --vh_col h_chain --vl_col l_chain -s