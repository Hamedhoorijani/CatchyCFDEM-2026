# movie of solvated peptide data

d = dump("files/dump.peptide")
d.unwrap()
p = pdbfile("files/peptide",d)

r = rasmol(p)
r.all()
