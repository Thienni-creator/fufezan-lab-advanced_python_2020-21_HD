import plotly.io as pio
pio.renderers.default = "browser"
from collections import deque
import requests

class Protein(object):
    def __init__(self, protein_id):
        """
        initializes an object of class protein for a given protein id

        Args:
            protein_id: Is the unique ID given to each protein

        """
        self.protein_id = protein_id

    def get_data(self):
        """
        downloads the fasta-file and extracts the 1-letter code of the amino acids
        for the protein object

        Returns:
            String containing the 1-letter code of each amino acid contained in the protein

        """
        url = 'https://www.uniprot.org/uniprot/' + self.protein_id + '.fasta?fil=reviewed:yes'
        r = requests.get(url)
        fasta_file_protein = 'C:/Users/ducth/PycharmProjects/thienni-creator_fufezan-lab-advanced_python_2020-21_HD' \
                             '/data/' + self.protein_id + '.fasta'
        with open(fasta_file_protein, 'wb') as file:
            file.write(r.content)
            file.close()

        with open(fasta_file_protein, 'r') as file:
            read_data = file.read()
            read_data_split = read_data.split('\n')
            sequences = ""
            for line in read_data_split:
                if not line.startswith('>'):
                    sequences += line
        return sequences

    def map(self, lookup:dict, aa_property:str, window_len = 1):
        """
        matches each amino acid in the sequence with its property and
        substitutes the values before averaging them for a sliding window

        Args:
            lookup: Dictionary containing
            aa_property: The property the protein has to be analyzed
            window_len: Number of amino acids which the property will be averaged over

        Returns:
            List containing the property for the amino acids after having been averaged
            for the given sliding window

        """
        mapping_dict = lookup[aa_property]
        seq = self.get_data()

        property_list = [mapping_dict[aa] for aa in seq]
        window = deque([], window_len)
        property_mean = []
        for aa_property in property_list:
            window.append(aa_property)
            window_mean = sum(window) / len(window)
            property_mean.append(window_mean)
            property_mean = [round(num, 2) for num in property_mean]
        return property_mean






