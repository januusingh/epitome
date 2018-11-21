# imports
from enum import Enum

######################################################
################### CONSTANTS ########################
######################################################
class Dataset(Enum):
    TRAIN = 1
    VALID = 2
    TEST = 3

    
# Regions in allpos.bed file that DeepSEA uses for training.
# Note that they are contiguous, so the total regions used 
# starts at _TRAIN_REGIONS[0] and ends at _TEST_REGIONS[1]+1
_TRAIN_REGIONS = [0, 2200000-1]
_VALID_REGIONS = [2200000, 2204000-1]
# 227512 for test (chr8 and chr9) for 227512 total rows
_TEST_REGIONS  = [2204000, 2204000 + 227512-1] # only up to chrs 8 and 9


    
# Modified regions in allpos.bed file that DeepSEA uses for training.
# We modified these because the validation set was too small.
_MODIFIED_TRAIN_REGIONS = np.r_[0:2200000,2400000:4200000]

# regions that should be taken from the train set for the validation set
_MODIFIED_VALID_FROM_TRAIN_REGIONS = np.r_[2200000:2400000,4200000:4400000]
# regions that should be taken from validation
_MODIFIED_VALID_FROM_VALID_REGIONS = np.r_[4200000:4208000]
# ORDER MATTERS!
_MODIFIED_VALID_REGIONS = np.concatenate([_MODIFIED_VALID_FROM_TRAIN_REGIONS, _MODIFIED_VALID_FROM_VALID_REGIONS])

# 227512 for test (chr8 and chr9) for 227512 total rows
_MODIFIED_TEST_REGIONS  = np.r_[2204000*2:2204000*2 + 227512*2] 



# DNase filepath dictionary
dnase_file_dict = {
  "A549": _ENCODE_DATA_PREFIX + "ENCFF001ARO_sorted_A549.bam",
  "HepG2": _ENCODE_DATA_PREFIX + "ENCFF224FMI_sorted_HepG2.bam",
  "K562": _ENCODE_DATA_PREFIX + "ENCFF271LGJ_sorted_K562.bam",
  'GM12878':_ENCODE_DATA_PREFIX + "ENCFF775ZJX_sorted_GM12878.bam",
  'H1-hESC':_ENCODE_DATA_PREFIX + "ENCFF571SSA_sorted_H1heSC.bam",
  'HeLa-S3':_ENCODE_DATA_PREFIX + "ENCFF783TMX_sorted_HeLaS3.bam",
  'HUVEC':_ENCODE_DATA_PREFIX + "ENCFF757PTA_sorted_HUVEC.bam",
  'GM12891':_ENCODE_DATA_PREFIX + "ENCFF070BAN_sorted_GM12891.bam",
  'MCF-7':_ENCODE_DATA_PREFIX + "ENCFF432OZA_sorted_MCF7.bam",
  'GM12892': _ENCODE_DATA_PREFIX + "ENCFF260LKE_sorted_GM12892.bam",
  'HCT-116': _ENCODE_DATA_PREFIX + "ENCFF291HHS_sorted_HCT116.bam"
}
