# pocket_hunting
Utilidades para evaluar pockets

# Pre-requisitos:
Es necesario bajarse CSA
https://www.ebi.ac.uk/thornton-srv/m-csa/api/homologues_residues.json

Hay que ir a UniProt y descargarse el archivo fasta del SwissProt humano
A ese archivo (nombrarlo uniprot_sprot_h.fasta) hay que correrle por línea de comando:
<pre><code>
makeblastdb -dbtype prot -in uniprot_sprot_h.fasta
</code></pre>

<pre><code>
#ejemplo :
python3 PocketHunt/matcheador_pdb_uniprot.py 5a46 A D 623
</code></pre>

<pre><code>
#ejemplo :
python3 PocketHunt/CSA_residuos.py P11362 D 623
</code></pre>

<pre><code>
#ejemplo :
python PocketHunt/active_site_detection.py 1br6.pdb PT1
</code></pre>


