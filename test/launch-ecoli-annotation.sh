#!/bin/bash

bacterial-annotator annotate -i data/Illumina_2x300bp_Ecoli_RayAssembly.fasta -g  U00096.3 -o Ecoli-RayAssembly-Annotation -f "$@"

