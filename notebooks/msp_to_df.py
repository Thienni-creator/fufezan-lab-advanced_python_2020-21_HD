import pandas as pd
from re import findall, match
from collections import namedtuple
import timeit


def parse_mod(mod_str):
    mod = namedtuple('Mod', [
        'x', 'y', 'z'
    ])

    x, y, z = mod_str.split(',')
    return mod(int(x), y, z)


def parse_header(header_line):
    spectrum = namedtuple('Spectrum', [
        'sequence',
        'charge',
        'mods',
        'collision_energy'
    ])

    mod_expression = r'\((\d+,[A-Za-z]+,[A-Za-z_\-]+)\)'

    header_expression = r'Name: (?P<sequence>[A-Z]+)\/' + \
        r'(?P<charge>\d)_' + \
        r'(?P<mods>\d(' + mod_expression + r')*)_' + \
        r'(?P<collisionenergy>\d+(\.\d+)?)eV'

    m = match(header_expression, header_line)

    sequence = m.group('sequence')
    charge = int(m.group('charge'))
    mods = [parse_mod(m) for m in findall(mod_expression, m.group('mods'))]
    collision_energy = float(m.group('collisionenergy'))

    return spectrum(sequence, charge, mods, collision_energy)


def parse_peak(line):
    mz, intensity, _ = line.split('\t')
    return float(mz), float(intensity)


def msp_to_df(input_file,
              max_seq_len=30,
              min_ce=36, max_ce=40,
              mz_min=135, mz_max=1400):
    """
    Function to read spectrum data from .msp file and convert to dataframe.
    Args:
        input_file (str): path to .msp file
        max_seq_len (int): maximum acceptable sequence length
        min_ce (int): minimum collision energy of spectra to be included in df
        max_ce (int): maximum collision energy of spectra to be included in df
        mz_min (int): lower boundary for m/z to be included in df
        mz_max (int): lower boundary for m/z to be included in df

    Returns:
        df (pd.DataFrame or np.array):   spectrum information within defined
                                         parameters [n_spectra, n_features]
        seqs (pd.DataFrame or np.array): sequences
    """
    def generate():
        bins = {}
        skip_spectrum = False
        for line in file:
            if line.startswith('Name:'):
                spectrum = parse_header(line)
                if min_ce <= spectrum.collision_energy <= max_ce and \
                   len(spectrum.sequence) <= max_seq_len:
                    bins['sequence'] = spectrum.sequence
                else:
                    skip_spectrum = True
            elif line.isspace():
                if not skip_spectrum:
                    yield bins.copy()
                bins.clear()
                skip_spectrum = False
            elif not skip_spectrum and match(r'\d+\.\d+\s\d+(\.\d+)?\s".*"', line):
                mz, intensity = parse_peak(line)
                if mz_min <= mz <= mz_max:
                    bins[round(mz)] = intensity

    with open(input_file, 'r') as file:
        df = pd.DataFrame(generate()).fillna(0.0)

    seqs = df.sequence
    df.drop(columns='sequence', inplace=True)

    return df.reindex(sorted(df.columns), axis=1), seqs


if __name__ == '__main__':
    file_name = '../data/cptac2_mouse_hcd_selected.msp'
    df, seqs = msp_to_df(file_name)
    print(df)
