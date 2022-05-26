#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio.PDB import *
import math
import pandas as pd
import sys
import os
from Geometry import *

__author__ = "Juan Manuel Prieto"

"""
    centroid_distance = 10 # distancia alrededor del centro de masa del ligando
"""

def H_Bonds(structure,ligand, centroid_distance , caso):
    Interaciones = {'acceptors':['GLU','ASP','ASN','GLN','SER','THR','TYR','HID','HIE'],
    'donors': ['ASN' ,'GLN' ,'ARG' , 'SER' ,'THR' , 'TYR' , 'TRP' , 'LYS' , 'HIE' , 'HID' , 'HIP'] }
    
    model = structure[0]
    for chain in model.get_chains():
        for residue in chain.get_residues():
            if residue.get_resname() == ligand:
                Ligando_Centro = list(center_of_mass(residue))

    Residuos_Interes = []
    atoms = []

    for residue in chain.get_residues():
        Residuo_Center = list(center_of_mass(residue))
        if (math.dist(Ligando_Centro, Residuo_Center)) < centroid_distance:
            if residue.get_resname() in Interaciones[caso]:
                Residuos_Interes.append([residue.get_resname(), residue.get_id()[1]])
                for atom in residue:
                    Res_name = residue.get_resname()
                    Res_id = residue.get_id()[1]
                    atom_name = atom.get_name()
                    Coor = list(atom.get_coord())
                    Serial = atom.get_serial_number()
                    atoms.append([Serial, Res_id, Res_name, atom_name, Coor, Residuo_Center])
                # sitio_activo.loc[len(sitio_activo.index)] = [Serial ,Res_id, Res_name, atom_name ,Coor,Residuo_Center ]
    return atoms


def active_site_residues(structure, ligand, centroid_distance):
    model = structure[0]

    for chain in model.get_chains():
        for residue in chain.get_residues():
            if residue.get_resname() == ligand:
                Ligando_Centro = list(center_of_mass(residue))

    Residuos_Interes = []
    atoms = []

    for residue in chain.get_residues():
        Residuo_Center = list(center_of_mass(residue))
        if (math.dist(Ligando_Centro, Residuo_Center)) < centroid_distance:
            Residuos_Interes.append([residue.get_resname(), residue.get_id()[1]])
            for atom in residue:
                Res_name = residue.get_resname()
                Res_id = residue.get_id()[1]
                atom_name = atom.get_name()
                Coor = list(atom.get_coord())
                Serial = atom.get_serial_number()
                atoms.append([Serial, Res_id, Res_name, atom_name, Coor, Residuo_Center])
                # sitio_activo.loc[len(sitio_activo.index)] = [Serial ,Res_id, Res_name, atom_name ,Coor,Residuo_Center ]
    return atoms


def as2table(active_site_residues,table_output):
    sitio_activo = pd.DataFrame(columns=['Serial', 'Pos', 'Residue', 'Atom', 'Coordenadas', 'Centro de Masa'])
    for atom in active_site_residues:
        sitio_activo.loc[len(sitio_activo.index)] = atom
    if table_output == False:
        pass
    else:
        sitio_activo.to_csv(table_output, index=False)
    return sitio_activo


def as2pdb(pdb_file, sitio_activo, pdb_out):
    pos = sitio_activo['Serial'].tolist()

    with open(pdb_out,"w") as h:
        Entrada = open(pdb_file, 'r').readlines()
        for lines in Entrada:
            if 'ATOM' in lines[0:4]:
                if int(lines[6:11]) in pos:
                    h.write(lines)
            elif 'HETATM' in lines[0:6]:
                if int(lines[6:11]) in pos:
                    h.write(lines)



if __name__ == "__main__":

    import argparse
    
    parser = argparse.ArgumentParser(description='Detects active site from PDB of a given ligand')
    parser.add_argument('pdb', action='store')
    parser.add_argument('ligand', action='store', help="PDB ligand code")
    parser.add_argument('-cd', '--centroid_distance', action='store', default=10, type=int)
    parser.add_argument("-to", '--table_output', action='store', type=str)
    parser.add_argument("-po",'--pdb_output', type=str, default="active_site.pdb")
    parser.add_argument("-bo",'--pdb_output_hbond', type=str, default=None)


    args = parser.parse_args()


    if not os.path.exists(args.pdb):
        sys.stderr.write(f"'{args.pdb}' does not exist")
        sys.exit(1)

    pdb_parser = PDBParser()

    try:
        structure = pdb_parser.get_structure("", args.pdb)
    except NameError:
        print('PDB Erroneo')
        exit()

    as_residues = active_site_residues(structure, args.ligand, args.centroid_distance)

    if args.table_output != None:
        sitio_activo = as2table(as_residues, args.table_output+'.csv')
    else:
        sitio_activo = as2table(as_residues, False)

    as2pdb(args.pdb,sitio_activo, args.pdb_output)


    if args.pdb_output_hbond != None:
        # Aceptores
        as_residues = H_Bonds(structure,args.ligand, args.centroid_distance , 'acceptors')
        sitio_activo_inter = as2table(as_residues, False)
        as2pdb(args.pdb,sitio_activo_inter, args.pdb_output_hbond+'_aceptor.pdb')
        # Dadores 
        as_residues = H_Bonds(structure,args.ligand, args.centroid_distance , 'donors')
        sitio_activo_inter = as2table(as_residues, False)
        as2pdb(args.pdb,sitio_activo_inter, args.pdb_output_hbond+'_donors.pdb')

    
    