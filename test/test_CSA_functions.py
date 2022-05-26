import sys

sys.path.insert(0, "/home/fginob1/Documents/Bitgenia/Doctorado/bioinfo_tesis/pocket_hunting")

from PocketHunt.CSA_functions import residuo_catalitico

csa = [{
'roles_summary': 'increase nucleophilicity, proton acceptor, proton donor, steric role',
'residue_sequences': [
           {"uniprot_id":"P11362","code":"Asp","resid":623}
	]
}]

result = residuo_catalitico("P11362", "D", 623,csa)

assert result and 'increase nucleophilicity, proton acceptor, proton donor, steric role' == result[0],"el aa deberia tener funcion catalitica y no se detecta"

result = residuo_catalitico("P11362", "E", 623,csa)

assert not result ,"no deberia pegar contra nada"

print("tudo bom")
