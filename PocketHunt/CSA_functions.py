import json
from Bio.SeqUtils import seq1

def residuo_catalitico(UniProtID, res, res_pos, csa):
    """
    Esta función toma como input un residuo y su localización en la secuencia canónica de UniProt
    y evalúa si el mismo figura en la base de datos del CSA como residuo catalítico
    
    Ejemplo:
    
    with open("catalytic_residues_homologues.json", "r") as handle:
        csa = json.load(handle)
    CSA_residuos("P11362", "D", 623,csa)
    """
    
    catalytic_functions = []
    
    # el CSA esta estructurado primeramente en función de roles catalíticos que pueden tener los residuos
    # catalíticos de una proteína, y dentro de esos roles describe qué secuencias hay documentadas, ya sea
    # por referencia o por homología
    for roles in csa:
        rol = roles['roles_summary']
        rol_seqs = roles['residue_sequences']
        for entry in rol_seqs:
            # chequeo la tríada UniProtID-residuo-posición
            if entry['uniprot_id'] == UniProtID and seq1(entry['code']) == res and entry['resid'] == res_pos:
                catalytic_functions.append(rol)

    return catalytic_functions

if __name__ == "__main__":
    
    # levanto la base de datos de CSA
    with open("catalytic_residues_homologues.json", "r") as handle:
        csa = json.load(handle)

    cs = residuo_catalitico("P11362", "D", 623, csa)

    print(cs)
