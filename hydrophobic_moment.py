#!/usr/bin/env python

"""
Calculates a set of properties from a protein sequence:
    - hydrophobicity (according to a particular scale)
    - mean hydrophobic dipole moment assuming it is an alpha-helix.
    - total charge (at pH 7.4)
    - amino acid composition
    - discimination factor according to Rob Keller (IJMS, 2011)

Essentially the same as HeliQuest (reproduces the same values).

Author:
  Joao Rodrigues
  j.p.g.l.m.rodrigues@gmail.com
"""

from __future__ import print_function

import argparse
import csv
import math
import os
import time

#
# Definitions
#
scales = {'Fauchere-Pliska': {'A':  0.31, 'R': -1.01, 'N': -0.60,
                              'D': -0.77, 'C':  1.54, 'Q': -0.22,
                              'E': -0.64, 'G':  0.00, 'H':  0.13,
                              'I':  1.80, 'L':  1.70, 'K': -0.99,
                              'M':  1.23, 'F':  1.79, 'P':  0.72,
                              'S': -0.04, 'T':  0.26, 'W':  2.25,
                              'Y':  0.96, 'V':  1.22},

          'Eisenberg': {'A':  0.25, 'R': -1.80, 'N': -0.64,
                        'D': -0.72, 'C':  0.04, 'Q': -0.69,
                        'E': -0.62, 'G':  0.16, 'H': -0.40,
                        'I':  0.73, 'L':  0.53, 'K': -1.10,
                        'M':  0.26, 'F':  0.61, 'P': -0.07,
                        'S': -0.26, 'T': -0.18, 'W':  0.37,
                        'Y':  0.02, 'V':  0.54},
          }
_supported_scales = list(scales.keys())

aa_charge = {'E': -1, 'D': -1, 'K': 1, 'R': 1}

#
# Functions
#
def assign_hydrophobicity(sequence, scale='Fauchere-Pliska'):  # noqa: E302
    """Assigns a hydrophobicity value to each amino acid in the sequence"""

    hscale = scales.get(scale, None)
    if not hscale:
        raise KeyError('{} is not a supported scale. '.format(scale))

    hvalues = []
    for aa in sequence:
        sc_hydrophobicity = hscale.get(aa, None)
        if sc_hydrophobicity is None:
            raise KeyError('Amino acid not defined in scale: {}'.format(aa))
        hvalues.append(sc_hydrophobicity)

    return hvalues


def calculate_moment(array, angle=100):
    """Calculates the hydrophobic dipole moment from an array of hydrophobicity
    values. Formula defined by Eisenberg, 1982 (Nature). Returns the average
    moment (normalized by sequence length)

    uH = sqrt(sum(Hi cos(i*d))**2 + sum(Hi sin(i*d))**2),
    where i is the amino acid index and d (delta) is an angular value in
    degrees (100 for alpha-helix, 180 for beta-sheet).
    """

    sum_cos, sum_sin = 0.0, 0.0
    for i, hv in enumerate(array):
        rad_inc = ((i*angle)*math.pi)/180.0
        sum_cos += hv * math.cos(rad_inc)
        sum_sin += hv * math.sin(rad_inc)
    return math.sqrt(sum_cos**2 + sum_sin**2) / len(array)


def calculate_charge(sequence, charge_dict=aa_charge):
    """Calculates the charge of the peptide sequence at pH 7.4
    """
    sc_charges = [charge_dict.get(aa, 0) for aa in sequence]
    return sum(sc_charges)


def calculate_discrimination(mean_uH, total_charge):
    """Returns a discrimination factor according to Rob Keller (IJMS, 2011)
    A sequence with d>0.68 can be considered a potential lipid-binding region.
    """
    d = 0.944*mean_uH + 0.33*total_charge
    return d


def calculate_composition(sequence):
    """Returns a dictionary with percentages per classes"""

    # Residue character table
    polar_aa = set(('S', 'T', 'N', 'H', 'Q', 'G'))
    speci_aa = set(('P', 'C'))
    apolar_aa = set(('A', 'L', 'V', 'I', 'M'))
    charged_aa = set(('E', 'D', 'K', 'R'))
    aromatic_aa = set(('W', 'Y', 'F'))

    n_p, n_s, n_a, n_ar, n_c = 0, 0, 0, 0, 0
    for aa in sequence:
        if aa in polar_aa:
            n_p += 1
        elif aa in speci_aa:
            n_s += 1
        elif aa in apolar_aa:
            n_a += 1
        elif aa in charged_aa:
            n_c += 1
        elif aa in aromatic_aa:
            n_ar += 1

    return {'polar': n_p, 'special': n_s,
            'apolar': n_a, 'charged': n_c, 'aromatic': n_ar}


def analyze_sequence(name=None, sequence=None, window=18, verbose=False):
    """Runs all the above on a sequence. Pretty prints the results"""

    if not sequence:
        raise Exception('Either I need glasses or there is no sequence.')

    if not name:
        name = 'Unnamed'

    w = window
    if w < 0:
        w = len(sequence)  # automatically set

    outdata = []  # for csv writing

    # Processing...
    seq_len = len(sequence)
    print('[+] Analysing sequence {} ({} aa.)'.format(name, seq_len))
    print('[+] Using a window of {} aa.'.format(w))
    for seq_range in range(0, seq_len):

        seq_w = sequence[seq_range:seq_range+w]
        if seq_range and len(seq_w) < w:
            break

        # Numerical values
        z = calculate_charge(seq_w)
        seq_h = assign_hydrophobicity(seq_w)
        av_h = sum(seq_h)/len(seq_h)
        av_uH = calculate_moment(seq_h)
        d = calculate_discrimination(av_uH, z)

        # AA composition
        aa_comp = calculate_composition(seq_w)
        n_tot_pol = aa_comp['polar'] + aa_comp['charged']
        n_tot_apol = aa_comp['apolar'] + aa_comp['aromatic'] + aa_comp['special']  # noqa: E501
        n_charged = aa_comp['charged']  # noqa: E501
        n_aromatic = aa_comp['aromatic']  # noqa: E501

        _t = [name, sequence, seq_range+1, w, seq_w, z, av_h, av_uH, d,
              n_tot_pol, n_tot_apol, n_charged, n_aromatic]
        outdata.append(_t)

        if verbose:
            print('  Window {}: {}-{}-{}'.format(seq_range+1, seq_range,
                                                 seq_w, seq_range+w))
            print('    z={:<3d} <H>={:4.3f} <uH>={:4.3f} D={:4.3f}'.format(z, av_h,  # noqa: E501
                                                                           av_uH, d))  # noqa: E501
            print('    Amino acid composition')
            print('      Polar    : {:3d} / {:3.2f}%'.format(n_tot_pol, n_tot_pol*100/w))  # noqa: E501
            print('      Non-Polar: {:3d} / {:3.2f}%'.format(n_tot_apol, n_tot_apol*100/w))  # noqa: E501
            print('      Charged  : {:3d} / {:3.2f}%'.format(n_charged, n_charged*100/w))  # noqa: E501
            print('      Aromatic : {:3d} / {:3.2f}%'.format(n_aromatic, n_aromatic*100/w))  # noqa: E501
            print()

    return outdata


def read_fasta_file(afile):
    """Parses a file with FASTA formatted sequences"""

    if not os.path.isfile(afile):
        raise IOError('File not found/readable: {}'.format(afile))

    sequences = []
    seq_name, cur_seq = None, None
    with open(afile) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith('>'):
                if cur_seq:
                    sequences.append((seq_name, ''.join(cur_seq)))
                seq_name = line[1:]
                cur_seq = []
            elif line:
                cur_seq.append(line)
    sequences.append((seq_name, ''.join(cur_seq)))  # last seq

    return sequences


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    i_opts = ap.add_mutually_exclusive_group(required=True)
    i_opts.add_argument('-s', '--sequence',
                        help='Sequence of amino acids to analyze')
    i_opts.add_argument('-f', '--seqfile',
                        help='File with sequences in FASTA format')
    ap.add_argument('-o', '--outfile',
                    help='File to write results in CSV format')
    ap.add_argument('-v', '--verbose', action='store_true',
                    help='Write information to screen as well')
    ap.add_argument('--scale', choices=_supported_scales,
                    help='Hydrophobicity scale to use')
    ap.add_argument('-w', '--window', default=18, type=int,
                    help=(
                        'AA window to use during analysis. Set to -1 to '
                        'automatically match the full-length of the sequence'))
    cmd = ap.parse_args()

    # File or Sequence?
    if cmd.sequence:
        all_data = analyze_sequence(sequence=cmd.sequence, window=cmd.window,
                                    verbose=cmd.verbose)
    elif cmd.seqfile:
        seq_list = read_fasta_file(cmd.seqfile)
        all_data = []
        for name, seq in seq_list:
            data = analyze_sequence(name=name, sequence=seq, window=cmd.window,
                                    verbose=cmd.verbose)
            all_data += data

    if not cmd.outfile:
        if cmd.seqfile:
            root, _ = os.path.splitext(cmd.seqfile)
        else:
            root = 'seq'
        outfn = root + '_hydrophobicity.txt'
    else:
        outfn = cmd.outfile

    print('[+] Writing results to {}'.format(outfn))
    if os.path.isfile(outfn):
        root, _ = os.path.splitext(outfn)
        outfn = root + '_' + time.strftime('%H%M%S%d%m%Y') + '.txt'
        print('  File already exists')
        print('  Writing to {}'.format(outfn))

    with open(outfn, 'w') as csvfile:
        _header = ['Name', 'Sequence', 'Window', 'Window Size', 'Sub-Sequence',
                   'Charge', 'Mean Hydrophobicity', 'Mean Hydrophobic Moment',
                   'Discrimination Factor', 'No. Polar AA', 'No. Apolar AA',
                   'No. Charged AA', 'No. Aromatic AA']

        writer = csv.writer(csvfile, dialect='excel')
        writer.writerow(_header)
        writer.writerows(all_data)
