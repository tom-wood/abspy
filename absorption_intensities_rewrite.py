fpath = "C:/Users/vgx18551/Documents/Data/"
fname = fpath + "absorption_intensities.txt"
new_fname = fpath + "macro_abs_intensities.txt"

full_text = []
with open(fname, 'r') as f:
    for line in f:
        l = line.strip()
        full_text.append(l)

new_text = "macro Cylinder_Absorption_Correction(murc, murv) {\n"
new_text += "#m_argu murc\n"
new_text += "If_Prm_Eqn_Rpt(murc, murv, min 0.1, max 20)\n"
mun = "prm #m_unique "
cev = "CeV(murc, murv)"
for l in full_text[1:]:
    if ";" in l:
        l += "\n"
    else:
        l += " "
    new_text += l.replace('prm !', mun).replace("mur", cev) 

new_text += "}"
with open(new_fname, 'w') as f:
    f.write(new_text)