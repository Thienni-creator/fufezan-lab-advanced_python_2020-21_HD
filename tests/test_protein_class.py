from classes import protein_class as pc
import pandas as pd


def test_get_data():
    protein1 = pc.Protein('P32249')
    gpcr = protein1.get_data()
    assert gpcr == "MDIQMANNFTPPSATPQGNDCDLYAHHSTARIVMPLHYSLVFIIGL" \
                   "VGNLLALVVIVQNRKKINSTTLYSTNLVISDILFTTALPTRIAYYA" \
                   "MGFDWRIGDALCRITALVFYINTYAGVNFMTCLSIDRFIAVVHPLR" \
                   "YNKIKRIEHAKGVCIFVWILVFAQTLPLLINPMSKQEAERITCMEY" \
                   "PNFEETKSLPWILLGACFIGYVLPLIIILICYSQICCKLFRTAKQN" \
                   "PLTEKSGVNKKALNTIILIIVVFVLCFTPYHVAIIQHMIKKLRFSN" \
                   "FLECSQRHSFQISLHFTVCLMNFNCCMDPFIYFFACKGYKRKVMRM" \
                   "LKRQVSVSISSAVKSAPEENSREMTETQMMIHSKSSNGK"


def test_map():
    property_dictionary = pd.read_csv("C:/Users/ducth/PycharmProjects/"
                                      "thienni-creator_fufezan-lab-advanced_python_2020-21_HD/data/"
                                      "amino_acid_properties.csv")
    property_dictionary = pd.DataFrame.to_dict(property_dictionary)
    one_letter_code = property_dictionary['1-letter code'].values()
    lookup = {}
    for pos, aa_property in enumerate(property_dictionary.keys()):
        if pos > 2:
            prop = property_dictionary[aa_property].values()
            prop_list = dict(zip(one_letter_code, prop))
            lookup[aa_property] = prop_list
    protein1 = pc.Protein('P32249')
    hydropathy_one = protein1.map(lookup, "hydropathy index (Kyte-Doolittle method)")
    assert hydropathy_one == [1.9, -3.5, 4.5, -3.5, 1.9, 1.8, -3.5, -3.5, 2.8, -0.7, -1.6, -1.6, -0.8,
                  1.8, -0.7, -1.6, -3.5, -0.4, -3.5, -3.5, 2.5, -3.5, 3.8, -1.3, 1.8, -3.2,
                  -3.2, -0.8, -0.7, 1.8, -4.5, 4.5, 4.2, 1.9, -1.6, 3.8, -3.2, -1.3, -0.8,
                  3.8, 4.2, 2.8, 4.5, 4.5, -0.4, 3.8, 4.2, -0.4, -3.5, 3.8, 3.8, 1.8, 3.8,
                  4.2, 4.2, 4.5, 4.2, -3.5, -3.5, -4.5, -3.9, -3.9, 4.5, -3.5, -0.8, -0.7,
                  -0.7, 3.8, -1.3, -0.8, -0.7, -3.5, 3.8, 4.2, 4.5, -0.8, -3.5, 4.5, 3.8,
                  2.8, -0.7, -0.7, 1.8, 3.8, -1.6, -0.7, -4.5, 4.5, 1.8, -1.3, -1.3, 1.8,
                  1.9, -0.4, 2.8, -3.5, -0.9, -4.5, 4.5, -0.4, -3.5, 1.8, 3.8, 2.5, -4.5,
                  4.5, -0.7, 1.8, 3.8, 4.2, 2.8, -1.3, 4.5, -3.5, -0.7, -1.3, 1.8, -0.4,4.2,
                  -3.5, 2.8, 1.9, -0.7, 2.5, 3.8, -0.8, 4.5, -3.5, -4.5, 2.8, 4.5, 1.8, 4.2,
                  4.2, -3.2, -1.6, 3.8, -4.5, -1.3, -3.5, -3.9, 4.5, -3.9, -4.5, 4.5, -3.5,
                  -3.2, 1.8, -3.9, -0.4, 4.2, 2.5, 4.5, 2.8, 4.2, -0.9, 4.5, 3.8, 4.2, 2.8,
                  1.8, -3.5, -0.7, 3.8, -1.6, 3.8, 3.8, 4.5, -3.5, -1.6, 1.9, -0.8, -3.9,
                  -3.5, -3.5, 1.8, -3.5, -4.5, 4.5, -0.7, 2.5, 1.9, -3.5, -1.3, -1.6, -3.5,
                  2.8, -3.5, -3.5, -0.7, -3.9, -0.8, 3.8, -1.6, -0.9, 4.5, 3.8, 3.8, -0.4,
                  1.8, 2.5, 2.8, 4.5, -0.4, -1.3, 4.2, 3.8, -1.6, 3.8, 4.5, 4.5, 4.5, 3.8,
                  4.5, 2.5, -1.3, -0.8, -3.5, 4.5, 2.5, 2.5, -3.9, 3.8, 2.8, -4.5, -0.7, 1.8,
                  -3.9, -3.5, -3.5, -1.6, 3.8, -0.7, -3.5, -3.9, -0.8, -0.4, 4.2, -3.5, -3.9,
                  -3.9, 1.8, 3.8, -3.5, -0.7, 4.5, 4.5, 3.8, 4.5, 4.5, 4.2, 4.2, 2.8, 4.2, 3.8,
                  2.5, 2.8, -0.7, -1.6, -1.3, -3.2, 4.2, 1.8, 4.5, 4.5, -3.5, -3.2, 1.9, 4.5,
                  -3.9, -3.9, 3.8, -4.5, 2.8, -0.8, -3.5, 2.8, 3.8, -3.5, 2.5, -0.8, -3.5, -4.5,
                  -3.2, -0.8, 2.8, -3.5, 4.5, -0.8, 3.8, -3.2, 2.8, -0.7, 4.2, 2.5, 3.8, 1.9,
                  -3.5, 2.8, -3.5, 2.5, 2.5, 1.9, -3.5, -1.6, 2.8, 4.5, -1.3, 2.8, 2.8, 1.8,
                  2.5, -3.9, -0.4, -1.3, -3.9, -4.5, -3.9, 4.2, 1.9, -4.5, 1.9, 3.8, -3.9,
                  -4.5, -3.5, 4.2, -0.8, 4.2, -0.8, 4.5, -0.8, -0.8, 1.8, 4.2, -3.9, -0.8,
                  1.8, -1.6, -3.5, -3.5, -3.5, -0.8, -4.5, -3.5, 1.9, -0.7, -3.5, -0.7, -3.5,
                  1.9, 1.9, 4.5, -3.2, -0.8, -3.9, -0.8, -0.8, -3.5, -0.4, -3.9]




