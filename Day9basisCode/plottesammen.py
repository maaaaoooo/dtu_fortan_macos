import numpy as np                                        # arrays
import matplotlib.pyplot as plt                           # plotting
from matplotlib.collections import PolyCollection         # tegner elementene som polygoner


def les_blokk(tekst, navn):                               # henter ut tallene i blokken "navn = [ ... ];"
    start = tekst.index(navn + " = [") + len(navn) + 4    # posisjonen rett etter "["
    slutt = tekst.index("];", start)                      # slutten på blokken
    return np.array([[float(v) for v in linje.split()]    # hver linje blir en rad med tall
                     for linje in tekst[start:slutt].strip().splitlines() if linje.strip()])

farger = {5: "gold", 30: "red", 95: "blue"}               # farge per x-verdi (x * 100)
xverdier = [5, 30, 95]                                    # x * 100, samme som i filnavnene
moder = [1, 2, 3]                                         # modene Fortran plottet
fig, akser = plt.subplots(3, 1, figsize=(3.0, 4.0)) 

for c, m in enumerate(moder):                             # én akse per modus
    ax = akser[c]                                         # aksen for denne moden
    xmin, xmax, ymin, ymax = 1e9, -1e9, 1e9, -1e9         # utstrekning, oppdateres under
    for xv in xverdier:                                   # tegn alle x-verdiene i samme akse
        fil = f"lbeam_x{xv:03d}_ploteig_data_{m:04d}.m"   # filnavnet Fortran lagde
        with open(fil) as f:                              # åpner filen og lukker den automatisk
            tekst = f.read()                              # leser hele filen
        X = les_blokk(tekst, "X")[:, :2]                  # nodekoordinater (x, y)
        IX = les_blokk(tekst, "IX").astype(int) - 1       # elementnoder, 0-indeksert
        D = les_blokk(tekst, "D").flatten()               # egenvektor (ux1, uy1, ux2, ...)
        D = D * np.sign(D[np.argmax(np.abs(D))])          # fikser fortegnet, ellers kan formen bli speilvendt
        Xdef = X + np.column_stack([D[0::2], D[1::2]])    # deformerte koordinater
        ax.add_collection(PolyCollection(Xdef[IX], facecolor=farger[xv], alpha=0.35,
                                         edgecolor=farger[xv], linewidth=0.3))   # gjennomsiktig polygon per element
        xmin, xmax = min(xmin, Xdef[:, 0].min()), max(xmax, Xdef[:, 0].max())   # utvider x-grensene
        ymin, ymax = min(ymin, Xdef[:, 1].min()), max(ymax, Xdef[:, 1].max())   # utvider y-grensene
    ax.set_xlim(xmin - 0.2, xmax + 0.2)                   # akseutstrekning i x, liten kant
    ax.set_ylim(ymin - 0.2, ymax + 0.2)                   # akseutstrekning i y, liten kant
    ax.set_aspect("equal")                                # lik skala på aksene
    ax.axis("off")                                        # fjerner akser og tall
    ax.set_title(f"Mode {m}", fontsize=8, pad=0)  

handles = [plt.Rectangle((0, 0), 1, 1, facecolor=farger[xv], edgecolor=farger[xv], alpha=0.5)
           for xv in xverdier]                            # én fargerute per x-verdi
labels = [f"x = {xv/100:.2f}" for xv in xverdier]         # kort tekst, tykkelsene forklares i bildeteksten
akser[0].legend(handles, labels, loc="lower right", fontsize=8, frameon=True)   # i den tomme delen av første rute
plt.tight_layout(h_pad=0.2) 
plt.savefig("modes_lbeam.png", bbox_inches="tight", pad_inches=0.05, dpi=300)   # lagrer figuren
plt.show()                                                # viser figuren                            # viser figuren