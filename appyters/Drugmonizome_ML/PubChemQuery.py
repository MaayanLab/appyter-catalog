import requests
import time
from ExponentialBackoff import ExponentialBackoff

PUBCHEM_BASE_URL = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound'

class PubChemQuery:
    """
    Class for making queries through the PubChem REST API.
    Uses exponential backoff to limit server request rates when 503 error is encountered.
    """
    
    backoff = ExponentialBackoff(min_value=0.2)

    @classmethod
    def make_query(cls, url):
        """
        Make a GET request from a specified URL with exponential backoff when a 503 error is encountered.
        The waiting time is halved upon each successful request and doubled upon each 503 response.

        Parameters:
            url (str): REST url to query
        
        Returns:
            (requests object): GET response
        """
        while True:
            time.sleep(cls.backoff.value())
            r = requests.get(url)
            throttling = r.headers['X-Throttling-Control']
            if ('Request Count status: Green' not in throttling) or ('Request Time status: Green' not in throttling) or ('too many requests' in throttling):
                cls.backoff.double()
                time.sleep(cls.backoff.value())
            if r.status_code != 503:
                cls.backoff.halve()
                break
            cls.backoff.double()
        return r

    @classmethod
    def query_by_name(cls, name, property):
        """
        Query PubChem API for drug with given name for a requested property.

        Parameters:
            name (str):     drug name
            property (str): property field as described by the PubChem API
        
        Returns:
            (requests object): GET response
        """
        return cls.make_query(f'{ PUBCHEM_BASE_URL }/name/{ name }/property/{ property }/TXT')

    @classmethod
    def query_by_smiles(cls, smiles, property):
        """
        Query PubChem API for drug with given SMILES for a requested property.

        Parameters:
            smiles (str):   drug SMILES
            property (str): property field as described by the PubChem API
        
        Returns:
            (requests object): GET response
        """
        return cls.make_query(f'{ PUBCHEM_BASE_URL }/smiles/{ smiles }/property/{ property }/TXT')

    @classmethod
    def name_to_inchi_keys(cls, name):
        """
        Query PubChem API for drug with given name for set of corresponding InChI Keys.

        Parameters:
            name (str): drug name
        
        Returns:
            (set): set of InChI Keys (str), empty if no matching drug is found
        
        Raises:
            RuntimeError: if the API response has status code other than 200 (successful) or 404 (no match)
        """
        r = cls.query_by_name(name, 'InChIKey')

        if r.status_code == 200:
            return set(r.content.decode('utf-8').split('\n')[:-1])
        elif r.status_code == 404:
            return set()
        else:
            raise RuntimeError(f'Unexpected HTTP status code { r.status_code } from PubChem response')

    @classmethod
    def name_to_smiles(cls, name, isomeric=True):
        """
        Query PubChem API for drug with given name for set of corresponding SMILES.

        Parameters:
            name (str): drug name
        
        Keyword argumetns:
            isomeric (bool): if True, returns isomeric SMILES
                             otherwise, returns canonical SMILES
        
        Returns:
            (set): set of SMILES (str), empty if no matching drug is found
        
        Raises:
            RuntimeError: if the API response has status code other than 200 (successful) or 404 (no match)
        """
        if isomeric:
            r = cls.query_by_name(name, 'isomericSMILES')
        else:
            r = cls.query_by_name(name, 'canonicalSMILES')

        if r.status_code == 200:
            return set(r.content.decode('utf-8').split('\n')[:-1])
        elif r.status_code == 404:
            return set()
        else:
            raise RuntimeError(f'Unexpected HTTP status code { r.status_code } from PubChem response')