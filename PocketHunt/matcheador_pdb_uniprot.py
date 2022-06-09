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
import argparse

def get_acc_sp(record):
    """
    Esta funcion permite que al indexar SwissProt las keys del diccionario sean los UniProtID
    """
    parts = record.split("|")
    assert len(parts) == 3 and parts[0] == "sp"
    return parts[1]

def matcheador_pdb_uniprot(cristal, cadena, res, res_pos, pdb_file_path, 
                           uniprot_db_path, uniprot_fastas, tmp="/tmp/"):
    """
    Esta función toma como input cristal, cadena, residuo y numeración del mismo en el cristal,
    además de la base de datos de fastas de UniProt tanto en formato archivo como indexada por BioPython
    y devuelve el UniProtID de la proteína, su residuo y su numeración base 1 en la secuencia de
    referencia de UniProt
    """

    # genero una estructura con freesasa que me permite iterar a través de los residuos con su numeración
    # propia del pdb (que puede no coincidir con la numeración de UniProt)
    structure_freesasa = freesasa.Structure(pdb_file_path)
    result = freesasa.calc(structure_freesasa)
    residueAreas = result.residueAreas()
    
    # voy a iterar a través de los residuos hasta llegar al residuo de mi interés (chequeo tanto por el
    # tipo de residuo como por su número), en el camino cuento cuántos residuos fui recorriendo,
    # el objetivo es saber cual es el número de aparición de mi residuo
    numero_residuo_cristal = 0
    for key in residueAreas[cadena]:
        numero_residuo_cristal += 1
        if residueAreas[cadena][key].residueType == one_to_three(res) and int(key) == res_pos:
            break
    
    # parseo el cristal 
    parser = PDBParser(PERMISSIVE=1)
    structure_id = cristal
    structure = parser.get_structure(structure_id, pdb_file_path)
    lista_residuos = []
    # lista donde voy a poner todos los residuos que provengan de la cadena donde está el match
    for model in structure:
        for chain in model:
            if chain.id == cadena:
                for residue in chain:
                    if residue.id[0] == " ":
                        lista_residuos.append(seq1(residue.get_resname()))
    
    # hago una cadena continua con todos los residuos que almacené y escribo un archivo temporal formato fasta
    # para hacer el blast contra UniProt
    seq_bruta_cristal = "".join(lista_residuos)
    with open(f"{tmp}/temp_file.fasta", "w") as out_handle:
        out_handle.write(">" + cristal + "\n" + seq_bruta_cristal)
        out_handle.close()
    
    # blast contra UniProt
    os.system(f"blastp -query {tmp}/temp_file.fasta -db {uniprot_db_path} -out {tmp}/pru_blast.xml -evalue 1e-5 -outfmt 5")
    
    # parseo la salida del blast contra UniProt, solo el primer hsp que debería ser el que más se parezca
    # a mi proteína de interés
    # me quedo con el UniProtID, la secuencia del query y el hit alineadas (con los gaps) y donde arranca
    # la coincidencia del query respecto de la secuencia de cristal que dimos como input
    blast_qresult = SearchIO.read(f"{tmp}/pru_blast.xml", "blast-xml")
    
    # interruptor que sirve solamente para cortar el proceso cuando se encuentre el residuo de interés
    interruptor = False
    # confío en que el mejor hit siempre va a ser la proteína misma
    # puede haber más de un hsp y puede que mi residuo de interés no caiga en el primero, por eso itero
    for hsp in blast_qresult[0]:
        if interruptor == True:
            break
        UniProtID = hsp.hit_id.split("|")[1]
        seq_aln_cristal = hsp.aln[0].seq
        seq_aln_uniprot = hsp.aln[1].seq
        start_query = hsp.query_range[0]
        # en primer lugar voy a recorrer la secuencia alineada del cristal y voy a ir contando gaps y residuos
        # hasta llegar al número de residuo que registré previamente
        numero_gaps_aln_cristal = 0
        numero_residuo_aln_cristal = start_query
        for char in seq_aln_cristal:
            if char == "-":
                numero_gaps_aln_cristal += 1
            else:
                numero_residuo_aln_cristal += 1
            if numero_residuo_aln_cristal == numero_residuo_cristal:
                # localización del residuo en la secuencia alineada, numeración base 0
                pos_aln_1 = numero_residuo_aln_cristal + numero_gaps_aln_cristal - start_query - 1
                # residuo en dicha localización pero sobre la secuencia de UniProt alineada
                residuo_uniprot_aln = seq_aln_uniprot[pos_aln_1]
                # chequeo que no sea un gap, sino no tiene sentido seguir
                if residuo_uniprot_aln != "-":
                    # voy a contar cuantos residuos hay sobre la secuencia alineada de UniProt en el
                    # alineamiento desde el comienzo hasta el residuo de interés (excluyo gaps)
                    numero_residuo_aln_uniprot = 0
                    for residuo_uniprot in seq_aln_uniprot[0:pos_aln_1 + 1]:
                        if residuo_uniprot != "-":
                            numero_residuo_aln_uniprot += 1
                    # me quedo mi secuencia de interés lista para meterla en un nuevo alineamiento
                    uniprot_seq_full = uniprot_fastas[UniProtID].seq
                    # alineamos la secuencia (parcial) que fue match contra el cristal, contra la
                    # secuencia completa que viene de UniProt
                    seq_1 = seq_aln_uniprot
                    seq_2 = uniprot_seq_full
                    aligner = Align.PairwiseAligner()
                    aligner.open_gap_score = -10 # muy penalizada la apertura para que me queden bloques grandes
                    aligner.extend_gap_score = -0.5
                    aligner.query_left_open_gap_score = 0
                    aligner.query_left_extend_gap_score = 0
                    aligner.query_right_open_gap_score = 0
                    aligner.query_right_extend_gap_score = 0
                    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
                    alignments = aligner.align(seq_1, seq_2)
                    alignment = alignments[0] #me quedo el primer alineamiento
                    # recojo las secuencias de los alineamientos con los gaps metidos en el medio
                    seq_uniprot_aln = str(alignment).split("\n")[0]
                    seq_full_uniprot_aln = str(alignment).split("\n")[2]

                    # cuento gaps y residuos sobre la secuencia parcial de UniProt hasta llegar
                    # al residuo de interés
                    numero_residuos_aln_fasta2 = 0
                    numero_gaps_aln_fasta2 = 0

                    for aa in seq_uniprot_aln:
                        if aa == "-":
                            numero_gaps_aln_fasta2 += 1
                        else:
                            numero_residuos_aln_fasta2 += 1
                        if numero_residuos_aln_fasta2 == numero_residuo_aln_uniprot:

                            # posición en base 0 del residuo en el alineamiento
                            pos_aln_2 = numero_residuos_aln_fasta2 + numero_gaps_aln_fasta2 - 1
                            # residuo sobre la secuencia completa de UniProt
                            residuo_seq_db_uniprot = seq_full_uniprot_aln[pos_aln_2]
                            interruptor = True
    
    # elimino los archivos temporales generados
    os.system(f"rm {tmp}/temp_file.fasta")
    os.system(f"rm {tmp}/pru_blast.xml")
    
    # retorno UniProtID, tipo de residuo y posición base 1
    return UniProtID, residuo_seq_db_uniprot, pos_aln_2+1

if __name__ == "__main__":
    
    warnings.filterwarnings("ignore")
    
    parser = argparse.ArgumentParser(description='Returns UniProtID, residue and its position (numeration 1-based) for the given residue from the crystal structure and chain')
    parser.add_argument('PDB', action='store', help = "PDBID, e.g. 5a46")
    parser.add_argument('Chain', action='store', help="Chain of the residue in the crystal structure e.g. A")
    parser.add_argument('Res', action='store', help="Residue of interest, one-letter code e.g. D")
    parser.add_argument('Res_pos', action='store', help="Position of the residue of interest, numeration 1-based e.g. 623", type=int)    
    parser.add_argument('-fastas', '--fasta_sequences', action='store', default="uniprot_sprot_h.fasta")

    args = parser.parse_args()
    
    # voy a buscar el archivo del cristal a la carpeta con todos los pdbs
    archivo = gzip.open("pdb" + args.PDB + ".ent.gz","rb")
    # leo el archivo y voy a escribir un archivo pdb para poder levantarlo luego con freesasa
    pdb = archivo.read()
    with open(args.PDB + ".pdb", "wb") as handle:
        handle.write(pdb)
        handle.close()

    # levanto la base de datos de secuencias de UniProt en formato fasta
    # las keys son los UniProtIDs, logrado gracias a la key_function
    uniprot_fastas = SeqIO.index_db(f"{args.fasta_sequences}.idx", args.fasta_sequences, "fasta", key_function=get_acc_sp)

    UniProtID, res, res_pos = matcheador_pdb_uniprot(args.PDB, args.Chain, args.Res, args.Res_pos,
                                                     args.PDB + ".pdb", args.fasta_sequences, uniprot_fastas)

    print(UniProtID, res, res_pos)
