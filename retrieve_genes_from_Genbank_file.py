#! /usr/bin/env python

####### Retrive a specified protein coding gene (CDS) from Genbank files and make fasta file for alignment ####
####### written by Filipe de Sousa on 02-02-2017, based on previous code by Cymon Cox  ####
####### usage: python retrive_Mt_genes.py <gene name 1> ####
####### run within a folder containing Genbank files ####



########################################################

import os
import sys
import copy
import Bio
#from gene_sampling import genes
#from taxon_sampling import pub_taxa
#from taxon_sampling import new_taxa
#from taxon_sampling import gene_accession_deletes
from BioSQL import BioSeqDatabase
#from Bio.Alphabet import IUPAC, Gapped
from Bio import SeqIO, GenBank
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, CompoundLocation


###################   ARGUMENTS    #####################

gene_name_1=sys.argv[1]



#######################################################



genbank_files = [x for x in  os.listdir(os.getcwd()) if ".gb" in x]   ### define where Genbank files are

sequences=[]  ### list to append gene sequences to



for genbankfile in genbank_files:
    print ("genbankfile is", genbankfile)
    record = SeqIO.read(genbankfile, "genbank")
    for feature in record.features:
        if feature.type == 'CDS'and "gene" in feature.qualifiers:
            if feature.qualifiers["gene"][0]== gene_name_1:
                print("found CDS", gene_name_1, "in", genbankfile, "!!!")
                sequence=feature.extract(record.seq)
                seq_r = SeqRecord(sequence, id=genbankfile, description="")
                sequences.append(seq_r)
            else:
                print("no genes here")

SeqIO.write(sequences, "%s_raw.fasta" % gene_name_1, "fasta")


print("DONE")
