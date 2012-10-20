#!/usr/bin/env python
"""
name: extractPeptidesFromFasta.py

usage: extractPeptidesFromFasta.py sequenceFile1 [sequenceFileN...]

commissioned by : Dr. Makoto Saito, 2012-10 .

authorship: adorsk, 2012-10, based on original perl script by Jeff Dusenberry.

description: This script extracts peptides from protein sequences stored in
FASTA files. For each file given as input it will generate a sorted 
list of rows containing: 
    <sequence_id>,<peptide_sequence>,<estimated_peptide_mass>
If multiple files are given, the results will be consolidated into a single file.
"""

"""
Imports and setup.
"""
import os, sys
import site
cur_dir = os.path.dirname(os.path.abspath( __file__ ))
site.addsitedir(os.path.join(cur_dir, "lib"))

def includePyteomics():
    """ 
    Some chicanery to include pyteomics. 
    Some parts of pyteomics require dependencies like 'matplotlib' or 'lxml', 
    which can be a PITA to install and distribute.
    We're not using those dependencies in this script, so we fake out
    pyteomics by pretending we have them.
    """
    class fakemodule(object):
        """ Mock stub for module imports. """
        @staticmethod
        def method(a, b):
            return a+b
    sys.modules["numpy"] = fakemodule
    sys.modules["lxml"] = fakemodule
    sys.modules["lxml"].etree = None
includePyteomics()
from pyteomics import fasta, parser, mass

"""
Process arguments.
"""
import argparse
argparser = argparse.ArgumentParser(description=('Extract peptides from FASTA'
                                              ' protein sequence files.'))
argparser.add_argument('fasta_files', nargs='+',
                    help=('List of FASTA files containing protein sequences.'
                          ' If no files are given, will use STDIN.'))
argparser.add_argument('--minmass', type=float, default=600.0,
                    help=('minimum peptide mass. Peptides with mass less'
                          ' than this will no be included in the output'))
argparser.add_argument('--maxmass', type=float, default=3500.0,
                    help=('maximum peptide mass. Peptides with mass greater'
                          ' than this will no be included in the output'))

"""
Main method.
"""
def main():
    args = argparser.parse_args()
    fasta_files = args.fasta_files

    # Only include peptides within these mass ranges.
    MIN_MASS = 600.0
    MAX_MASS = 3500.0

    # Read sequence records from files.
    sequence_records = []
    for fasta_file in fasta_files:
        sequence_records.extend(read_fasta_sequences(fasta_file))

    # Make a peptide info database.
    peptides_db = {}
    for sequence_record in sequence_records:
        for peptide in sequence_record['peptides']:
            if peptide not in peptides_db:
                peptides_db[peptide] = get_peptide_data(peptide)

    # Assemble data records.
    output_records = []
    for sequence_record in sequence_records:
        for peptide in sequence_record['peptides']:
            peptide_data = peptides_db[peptide]
            if MAX_MASS > peptide_data['mass'] > MIN_MASS:
                output_records.append({
                    'sequence_id': sequence_record['id'],
                    'peptide': peptide,
                    'peptide_mass': peptide_data['mass']
                })

    # Output the data.
    output_fields = ['sequence_id', 'peptide', 'peptide_mass']
    for output_record in output_records:
        output_row = [str(output_record[field]) for field in output_fields]
        print ','.join(output_row)

"""
Helper methods.
"""

def get_peptide_data(peptide):
    """ Get data for a given peptide. """
    peptide_data = {'sequence': peptide}
    peptide_data['parsed_sequence'] = parser.parse(
        peptide,
        show_unmodified_termini=True # keep the termini, for mass calculations.
    )
    peptide_data['mass'] = mass.calculate_mass(
        peptide_data['parsed_sequence']
    )
    return peptide_data


def read_fasta_sequences(fasta_file):
    """ Read sequence records from a FASTA file. """
    sequence_records = []
    for description, sequence in fasta.read(fasta_file):
        # Initialize sequence record with sequence string.
        sequence_record = {'sequence': sequence}

        # Get sequence info.
        description_parts = description.split()
        sequence_record['id'] = description_parts[0]

        # Get the sequence's peptides.
        sequence_record['peptides'] = parser.cleave(
            sequence, 
            parser.expasy_rules['trypsin'],
            missed_cleavages=1 #max no. of missed cleavages.
        )

        # Save the sequence record, keyed by the id.
        sequence_records.append(sequence_record)

    return sequence_records

if __name__ == '__main__':
    main()
