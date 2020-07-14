import io
import requests
import os
import pandas as pd

BASE_URL = 'https://appyters.maayanlab.cloud/storage/Drugmonizome_ML/SEP-L1000/'
DATASETS = ['GO_transformed_signatures_PAEA.csv.gz',
            'meta_SMILES_InChIKey_drugmonizome.csv.gz',
            'LINCS_Gene_Experssion_signatures_CD.csv.gz',
            'MLPCN_morplological_profiles.csv.gz',
            'MACCS_bitmatrix.csv.gz']

class SEPL1000(object):
    """Class for retrieving SEP-L1000 data stored in an online repository.
    """

    @classmethod
    def download(cls, datasets=None):
        """Download datasets from the SEP-L1000 data repository and store locally in directory 'SEPL1000'.
        """

        if datasets is None:
            datasets = DATASETS

        for filename in datasets:
            if not os.path.exists('SEPL1000'):
                os.mkdir('SEPL1000')
            filepath = 'SEPL1000/' + filename

            if not os.path.isfile(filepath):

                url = BASE_URL + filename
                response = requests.get(url, stream=True)

                if response.status_code != 200:
                    raise Exception('This should not happen. Filename may be invalid.')

                cls._download_file(response, filepath)

            yield filepath

    @classmethod
    def download_df(cls, datasets=None, **kwargs):
        for file in cls.download(datasets):
            yield pd.read_csv(file, **kwargs)

    @classmethod
    def _download_file(cls, response, filepath):
        """Saves the content of a response request to a specified file.
        """

        with open(filepath, 'wb') as outfile:
            for chunk in response.iter_content(chunk_size=1024):
                outfile.write(chunk)