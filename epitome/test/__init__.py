import os
import sys
import tempfile
import unittest
from epitome.models import *
from epitome import GET_DATA_PATH

class EpitomeTestCase(unittest.TestCase):

	def getValidData(self):
		data_path = GET_DATA_PATH()
		x = os.path.join(data_path, 'valid.npz')
		return scipy.sparse.load_npz(x).toarray()

	def getFeatureData(self, eligible_assays, eligible_cells):
		# returns matrix, cellmap, assaymap
		return get_assays_from_feature_file(eligible_assays = eligible_assays,
				eligible_cells = eligible_cells, min_cells_per_assay = 3, min_assays_per_cell = 1)

	def makeSmallModel(self):

		sparse_matrix = self.getValidData()
		data = {Dataset.TRAIN: sparse_matrix, Dataset.VALID: sparse_matrix, Dataset.TEST: sparse_matrix}


		eligible_cells = ['K562','HepG2','H1','A549','HeLa-S3']
		eligible_assays = ['DNase','CTCF']
		matrix, cellmap, assaymap = self.getFeatureData(eligible_assays, eligible_cells)

		return VLP(list(eligible_assays),
			test_celltypes = ['K562'],
			matrix = matrix,
			assaymap = assaymap,
			cellmap = cellmap,
			data = data)


	def tmpFile(self):

		tempFile = tempfile.NamedTemporaryFile(delete=True)
		tempFile.close()
		return tempFile.name