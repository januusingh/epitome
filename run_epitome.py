################################################################
############ Runs Epitome on Single Transcription Factor #####################
################################################################


import argparse
import os
from epitome.models import *
from epitome.functions import *
from epitome.viz  import *

from epitome.constants import *
import yaml
import subprocess
from timeit import default_timer as timer
# Parser for user specific locations
parser = argparse.ArgumentParser(description='Run Epitome on single TF.')

parser.add_argument('--epitome_data_path', help='path to epitome data')
parser.add_argument('--results_path', help='Path to save models and results to')

# results_path = parser.parse_args().results_path
# epitome_data_path = parser.parse_args().epitome_data_path

results_path = "results"
epitome_data_path = "data/epitome_data" 
feature_path = os.path.join(epitome_data_path, "feature_name")
#motif_data = "epitome/GATA2_mat.npz"
#motif_inds = "epitome/tf_ind.csv"

# if parser.parse_args().bam_data_path is not None:
#     bam_data_path = parser.parse_args().bam_data_path
# else:
#     bam_data_path = defcom_data_path

# create user directories if they do not exist
epitome_results_dir = os.path.join(results_path, "epitome_results")
model_dir = os.path.join(results_path, "epitome_models")
# tf_model_dir = os.path.join(results_path, "epitome_GATA2_models")

if not os.path.exists(epitome_results_dir):
    os.makedirs(epitome_results_dir)
if not os.path.exists(model_dir):
    os.makedirs(model_dir)
# if not os.path.exists(tf_model_dir):
#     os.makedirs(tf_model_dir)
# load in configuration

# load in data for epitome
#train_data, valid_data, test_data = load_epitome_data(epitome_data_path)#config['epitome_data_dir'])
train_data = scipy.sparse.load_npz(os.path.join(epitome_data_path, 'train.npz')).toarray()
valid_data = scipy.sparse.load_npz(os.path.join(epitome_data_path, 'valid.npz')).toarray()
test_data = scipy.sparse.load_npz(os.path.join(epitome_data_path, 'test.npz')).toarray()
data = {Dataset.TRAIN: train_data, Dataset.VALID: valid_data, Dataset.TEST: test_data}
# # load motif data for epitome TF model
# motifmat = np.load(motif_data)
# motifmap = pd.read_csv(motif_inds, header=None)

# Get TFs to train on
#TFs = list(map(lambda x: x[0], EPITOME_AND_DEFCOM_TFS))
TFs = 'CTCF'

matrix, cellmap, assaymap = get_assays_from_feature_file(feature_path,
                                                         eligible_assays = None,
                                  eligible_cells = None, min_cells_per_assay = 2, min_assays_per_cell= 2) #10)

# plot_assay_heatmap(matrix, cellmap, assaymap)
print("Generated non-sparse data matrix with %s cell lines" % len(cellmap))
print("Cellmap shape ", list(cellmap))

all_data = np.concatenate((data[Dataset.TRAIN], data[Dataset.VALID], data[Dataset.TEST]), axis=1)
eval_results_df = pd.DataFrame(columns=['transcription_factor', 'query_cell', 'auROC', 'auPRC'])
# iterate through 4 cell lines and train a model holding out each.
# then evaluate each of the 4 models on the single TF
query_cell = 'K562' #'T47D'
# Original Model
model = VLP(['CTCF', 'SMC3', 'RAD21'])
model.train(5000) # train for 5000 iterations

#epitome_model = VLP(data,
#               [query_cell],
#               matrix,
#               assaymap,
#               cellmap,
#               None,  #motifmat,
#               None, #motifmap,
#               shuffle_size=2,
#               batch_size=64)
#start = timer()
#epitome_model.train(2000)
#end = timer()
#print('epitome train: %f' % (end - start))
model_path = os.path.join(model_dir, query_cell)
model.save(model_path)
start = timer()
model_results = model.test(10000, calculate_metrics=True)
end = timer()
print('epitome test %f' % (end-start))
print('Model auROC: %s. Model auPRC: %s.' % (model_results['auROC'], model_results['auPRC'])) 
eval_results_df = eval_results_df.append({ 
   'transcription_factor' : 'CTCF_SMC3_RAD21',
   'query_cell' : query_cell,
   'auROC' : model_results['auROC'],
   'auPRC' : model_results['auPRC']}, 
    ignore_index=True)
eval_results_df.to_csv(os.path.join(epitome_results_dir,'epitome_original_epitome_data.csv'), sep="\t")


# Model with TF data
#TFs = 'CTCF'
#assays = [TFs] + ['DNase']
#matrix, cellmap, assaymap = get_assays_from_feature_file(feature_path,
#                                                         eligible_assays = TFs,
#                                  eligible_cells = None, min_cells_per_assay = 2, min_assays_per_cell= 2) #10)

#epitome_tf_model = VLP(data,
#               [query_cell],
#               matrix,
 #              assaymap,
 #              cellmap,
  #             motifmat,
  #             motifmap,
#               shuffle_size=2,
#               batch_size=64)
#start = timer()
#epitome_tf_model.train(2000)
#end = timer()
#print('epitome %s train: %f' % (TFs, end - start))
#model_path = os.path.join(tf_model_dir, query_cell)
#epitome_tf_model.save(model_path)
#start = timer()
#model_results = epitome_tf_model.test(10000, calculate_metrics=True)
#end = timer()
#print('epitome %s test %f' % (TFs, end-start))
#print('Model auROC: %s. Model auPRC: %s.' % (model_results['auROC'], model_results['auPRC']))
#eval_results_df = eval_results_df.append({
#   'transcription_factor' : TFs,
#   'query_cell' : query_cell,
#   'auROC' : model_results['auROC'],
#   'auPRC' : model_results['auPRC']},
#    ignore_index=True)
#eval_results_df.to_csv(os.path.join(epitome_results_dir,'epitome_CTCF_epitome_data.csv'), sep="\t")
# for query_cell_EPITOME_NAME, query_cell_DEFCOM_NAME in EPITOME_AND_DEFCOM_CELLS:
#
#    cell_model_dir = os.path.join(model_dir, "query_" + query_cell_EPITOME_NAME)
#    if not os.path.exists(cell_model_dir):
#        os.makedirs(cell_model_dir)
#
#    # Make sure checkpoint index file exists
#    if os.path.isfile(os.path.join(cell_model_dir,"saved_model.pb")):
#        print("Found saved model at %s" % cell_model_dir)
##        epitome_model = VLP(data = data, checkpoint=cell_model_dir)
#
#    else:
#        print("Invalid model path %s. Training new model..." % cell_model_dir)
#        # Train model from scratch
#        epitome_model= VLP(data,
#                    [query_cell_EPITOME_NAME], # cell line reserved for testing
#                    matrix,
#                    assaymap,
#                    cellmap,
#                    motifmat,
#                    motifmap, 
#                    shuffle_size=2,
#                    batch_size=64)
#        epitome_model.train(5000)
#
#        # Save model
#        epitome_model.save(cell_model_dir)
#    # Test model
#    model_results = epitome_model.test(100)
#    print('Model auROC: %s. Model auPRC: %s.' % (model_results['auROC'], model_results['auPRC']))
#    eval_results_df = eval_results_df.append({
#        'transcription_factor' : TF,
#	'query_cell' : query_cell_EPITOME_NAME, 
#	'auROC' : model_results['auROC'],
#	'auRPC' : model_results['auPRC']}, ignore_index=True)
#
# eval_results_df.to_csv(os.path.join(epitome_results_dir,'epitome_allresults_defcomData.csv'), index=False)


#eval_results_df = pd.DataFrame(columns=['model','training_cell','transcription_factor','query_cell','auROC','auPR'])

# for all the query cell types
#for query_cell_epitome_name, query_cell in EPITOME_AND_DEFCOM_CELLS:
    # for all the transcription factors
 #   for tf_epitome_name, tf in EPITOME_AND_DEFCOM_TFS:
        # get the files we'll need
 #       prediction_results_file = os.path.join(epitome_results_dir, '{}_{}_{}.csv'.format('joint', tf, query_cell))
 #       pos_file = os.path.join(defcom_data_path,'{}_{}_pos_test.bed'.format(tf, query_cell))
        # try just incase we are missing data for one
 #       try:
            # evaluate the model
 #           auROC, auPR = evaluateEpitomeResults(tf_epitome_name, prediction_results_file, pos_file)
            # append the results
 #           eval_results_df = eval_results_df.append({'model':'epitome',
 #                       'training_cell': 'joint',
 #                       'transcription_factor': tf,
 #                       'query_cell': query_cell,
 #                       'auROC': auROC,
 #                       'auPR': auPR}, ignore_index=True)
        # if something goes wrong we should know about it
 #       except Exception as e:
 #           print('Unable to score predictions for TF: {} and query cell type {}. Exception: {}'.format(tf, query_cell, e))

#eval_results_df.to_csv(os.path.join(epitome_results_dir,'epitome_allresults_defcomData.csv'), index=False)
