#!/bin/bash

bacterial-annotator annotate -i data/Illumina_2x300bp_Ecoli_RayAssembly.fasta\
		    -g U00096.3\
		    --remotedb swissprot\
		    -o Ecoli-RayAssembly-Annotation-REF+NCBIremote\
		    -f\
		    --pidentity 60\
		    --minlength 600\
		    --meta
