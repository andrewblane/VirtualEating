from __future__  import print_function
import os
import Bio
import sys
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio import Restriction
from Bio.Restriction import *
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqFeature
from Bio.SeqFeature import *
import subprocess

def generate_score_db(genome_directory = "", database_name= ""):
    ###
    '''
    Takes a directory containing genomic data in FASTA format and converts it
    to a BLAST database in format expected by EATING functions.

    Note that genomic data is often downloaded in a zipped format; unzip all
    files before use (using, for example, the command 'gzip -d *.gz' from within
    the directory containing the zipped files)
    '''
    #Directory with your files
    if genome_directory == "":
        # OS
        genome_directory = raw_input("Enter directory containing genome FASTA files:")

    if database_name == "":
        database_name = raw_input("Pick a short name for your BLAST database (e.g. hg19):")

    # Concatenate individual FA files
    file_list = subprocess.Popen(str("ls " + genome_directory), stdout=subprocess.PIPE, shell = True)
    if "_all.fa" in file_list.communicate()[0]:
        print("Found concatenated file ending in \"_all.fa\"")
    else:
        os.system(str("cat " + genome_directory + "/*.fa > " + genome_directory + "/" + database_name + "_all.fa"))

    '''
    Example:
    makeblastdb -in LAEVIS_7.1.repeatmasked.fa -dbtype nucl -parse_seqids -out hg19 -title 'hg19'
    /Volumes/2TB\ HFS/MBP\ 750GB/andypy/Genomic\ Data/hg19/
    -parse_seqids
    '''

    genome_directory_nonescaped = []
    for char in genome_directory:
        if char == "\\":
            None
        else:
            genome_directory_nonescaped.append(char)
    genome_directory_nonescaped = "".join(genome_directory_nonescaped)

    genome_path_b = str("\"\\\""+genome_directory_nonescaped+ database_name + "_all.fa\\\"\"")
    blast_db_generation_command = str("makeblastdb -in " + genome_path_b + " -dbtype nucl -out " + database_name + " -title " + database_name + " -parse_seqids")
    print(blast_db_generation_command)
    blast_command_output = subprocess.Popen(blast_db_generation_command, stdout=subprocess.PIPE, shell = True)
    print(blast_command_output.communicate())
