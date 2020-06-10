"""
Example script to generate an aggregated motif matrix for CTCF, REST, ATF3.
"""

import os
import pandas as pd
from epitome.motif_functions import *

queried_tfs = ['CTCF', 'REST', 'ATF3']
motif_save_dir = "data/motif_data"
suffix = "ctcf_rest_aft3_hoco"

generate_motif_matrix(queried_tfs, 
                     "/data/yosef2/DC/footprintPrediction/geneInductionPipeline/TF_motifs/HOCOMOCO/pwmscanOut",
                     "data/epitome_data/all.pos.bed",
                     motif_save_dir, 
                     suffix)

# Agg. motif sites for all data per TF
motifmat, motifmap = unique_motif_matrix(os.path.join(motif_save_dir, suffix + "_motifmat.npz"), 
                                         os.path.join(motif_save_dir, suffix + "_motifmap.csv"),
                                         motif_save_dir)
