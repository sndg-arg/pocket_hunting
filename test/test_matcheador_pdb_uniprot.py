import sys
from Bio import SearchIO
from Bio import AlignIO
import freesasa
from Bio.PDB.PDBParser import PDBParser
from Bio import Align
from Bio.Align import substitution_matrices
import os
from Bio import SeqIO
import gzip
import json
from Bio.SeqUtils import seq1
from Bio.PDB.Polypeptide import one_to_three
import warnings

warnings.filterwarnings("ignore")

sys.path.insert(0, "/home/fginob1/Documents/Bitgenia/Doctorado/bioinfo_tesis/pocket_hunting")

from PocketHunt.matcheador_pdb_uniprot import matcheador_pdb_uniprot
from PocketHunt.matcheador_pdb_uniprot import get_acc_sp

prueba_fasta = SeqIO.index_db("testfasta.fasta.idx", "testfasta.fasta", "fasta", key_function=get_acc_sp)

UniProtID, res, res_pos = matcheador_pdb_uniprot("5a46", "A", "D", 623, "5a46.pdb", "testfasta.fasta", prueba_fasta)

assert UniProtID == "P11362" and res == "D" and res_pos == 623, "el aa deberia corresponder a D623 de P11362"

print("tudo bom")

UniProtID, res, res_pos = matcheador_pdb_uniprot("5a46", "A", "H", 621, "5a46.pdb", "testfasta.fasta", prueba_fasta)

assert res != "D" and res_pos != 623, "no deberia pegar contra nada"

print("tudo bom")
