r"""

Functions for motif data.

================
Helper functions
================
  generate_motif_matrix
  unique_motif_matrix
"""

# Imports
import os
import sys
import pandas as pd
import numpy as np
import pybedtools
import argparse

def generate_motif_matrix(queried_tfs, 
                          motif_data_dir, 
                          all_regions_file, 
                          motif_save_dir, 
                          motif_save_suf):
    """
    Generates a (n x m) matrix, where n is the number of files examined 
    and m is the number of sites in the all_regions_file, and saves it to
    the motif_save_dir as motifmat. Also generates and saves a corresponding 
    motifmap which maps each row in motifmat to its TF name. 

    :param tfs: list of transcription factors.
    :param motif_data_dir: str directory where the motif data is stored. 
    :param all_regions_file: bed file containing a row for each genomic region.
    :param motif_save_dir: str directory to save motif matrix and map.
    :param motif_save_suf: str suffix to add to file names when saving matrix and map.

    :return
    """
    assert isinstance(queried_tfs, list), "queried_tfs must be a list."
    
    epitome_tfs_file = 'data/epitome_data/feature_name'
   
    # Output file names
    motifmat_file = os.path.join(motif_save_dir, motif_save_suf + "_motifmat")
    motifmap_file = os.path.join(motif_save_dir, motif_save_suf + "_motifmap.csv")
    missing_tfs_file = os.path.join(motif_save_dir, 'missing_tfs.csv')
    print('Output for motifmat file at %s ...' % (motifmat_file))
    print('Output for motifmap file at %s ...' % (motifmap_file))

    # Unique list of TFs in Epitome Dataset
    epitome_tfs = pd.read_table(epitome_tfs_file, header=None, skiprows=[0], sep='|')[1].unique()
    
    # List of TFs in motif_data_dir not in EPITOME_DATA
    missing_tfs, motifmap = [], []
    
    # Comparison Bed File
    a = pybedtools.BedTool(all_regions_file)
    num_motifs = 3268840 # a.count()
    motifmat = np.empty((0,num_motifs))
    has_found = False
    
    # Iterate through each motif file in motif_data_dir and check if the 
    # file's TF is the same as the queried TF
    for motif_f in os.listdir(motif_data_dir):
        motif_path = os.path.join(motif_data_dir, motif_f)
        file_tf = motif_f.split('_')[0]

        # Check if current TF in list of epitome TFs
        if (file_tf not in epitome_tfs):
            missing_tfs.append(file_tf)
            pd.Series(missing_tfs).to_csv(missing_tfs_file, header=False)
            continue
            
        elif os.path.isfile(motif_path) and (file_tf in queried_tfs):
            # For HOCOMOCO dataset only look at the sorted bed files.
            if ("HOCOMOCO" in motif_path) and ("sorted" not in motif_path):
                continue
                
            print('Looking at %s ...' % (motif_path))
            b = pybedtools.BedTool(motif_path)

            # Left intersect on all epitome binding sites and current file binding sites
            # Note: "c=True" gives the counts for each time the site is intersected in the file
            intr = a.intersect(b, c=True)
            overlap_counts = intr.to_dataframe()['name']

            # Turn counts into a 0 or 1 binary value of whether it intersected
            overlap_counts[overlap_counts > 0] = 1
            
            # Assert dimensions of TF analysis
            assert(overlap_counts.size == num_motifs)

            motifmat = np.append(motifmat, [overlap_counts], axis=0)
            motifmap.append(file_tf)

            # Saves data motif matrix and map at each iteration
            np.savez_compressed(motifmat_file, tf=motifmat)
            pd.Series(motifmap).to_csv(motifmap_file, header=False)

            # Indicate that at least 1 motif file has been found
            has_found = True

    # Check if atleast 1 motif file has matched with given motif
    if (not has_found):
        print("No motif sites found for transcription factor..")
        
    return

def unique_motif_matrix(motifmat_file, motifmap_file, save_dir=None):
    """
    Generates a (j x m) matrix, where j is the number of unique TFs
    and m is the number of sites in the all_regions_file. Also generates 
    and saves a motifmap which maps each row in motifmat to its TF name. 

    :param motifmat_file: str filename where the motifmat is stored.
    :param motifmap_file: str filename where the motifmap is stored. 
    :param save: str directory to store the matrix/map. Won't save if None.

    :return tuple (summed motifmat, summed motifmap)
    """
    # Load in data
    motifmat = np.load(motifmat_file)["tf"]
    motifmap = pd.read_csv(motifmap_file, header=None).rename(columns={0:"Index", 1:"TF"})
    
    # Unique TF's in motifmap
    unique_tf = motifmap["TF"].unique()
    print("Number of unique TFs: %s" % (unique_tf.shape))

    motifmat_sum = np.empty((0, 3268840))
    motifmap_sum = []

    for tf in unique_tf:
        tf_index = motifmap[motifmap["TF"] == tf]["Index"]
        tf_motif = motifmat[tf_index,:]
        tf_motif_sum = np.sum(tf_motif, axis=0)
        tf_motif_sum[tf_motif_sum > 0] = 1
        motifmat_sum = np.append(motifmat_sum, [tf_motif_sum], axis=0)
        motifmap_sum = np.append(motifmap_sum, tf)
    
    motifmap_sum = pd.Series(motifmap_sum)
        
    print("Summed motifmat shape: %s" % (str(motifmat_sum.shape)))
    print("Summed motifmap shape: %s" % (str(motifmap_sum.shape)))
    
    if save_dir is not None:
        motifmat_sum_file = os.path.join(save_dir, "unique_motifmat")
        motifmap_sum_file = os.path.join(save_dir, "unique_motifmap.csv")
        
        np.savez_compressed(motifmat_sum_file, tf=motifmat_sum)
        motifmap_sum.to_csv(motifmap_sum_file, header=False)
        
        print('Output for motifmat file at %s ...' % (motifmat_sum_file))
        print('Output for motifmap file at %s ...' % (motifmap_sum_file))
    
    return motifmat_sum, motifmap_sum
    
