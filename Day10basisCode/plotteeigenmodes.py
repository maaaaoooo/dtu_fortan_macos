import numpy as np                                        # arrays
import matplotlib.pyplot as plt                           # plotting
from matplotlib.collections import PolyCollection         # tegner elementene som polygoner

def les_blokk(tekst, navn):                               # henter ut tallene i blokken "navn = [ ... ];"
    start = tekst.index(navn + " = [") + len(navn) + 4    # posisjonen rett etter "[" 
    slutt = tekst.index("];", start)                      # slutten på blokken
    return np.array([[float(v) for v in linje.split()]    # hver linje blir en rad med tall
                     for linje in tekst[start:slutt].strip().splitlines() if linje.strip()])

xverdier = [5, 30, 95]                                   # x * 100, samme som i filnavnene
moder = [1, 2, 3]                                         # modene Fortran plottet
fig, akser = plt.subplots(3, 3, figsize=(12, 9))          # 3 x-verdier (rader) og 3 moder (kolonner)

for r, xv in enumerate(xverdier):                         # rad = x-verdi
    for c, m in enumerate(moder):                         # kolonne = modus
        fil = f"lbeam_x{xv:03d}_ploteig_data_{m:04d}.m"   # filnavnet Fortran lagde
        tekst = open(fil).read()                          # leser hele filen
        X = les_blokk(tekst, "X")[:, :2]                  # nodekoordinater (x, y)
        IX = les_blokk(tekst, "IX").astype(int) - 1       # elementnoder, 0-indeksert
        D = les_blokk(tekst, "D").flatten()               # egenvektor (ux1, uy1, ux2, ...)
        Xdef = X + np.column_stack([D[0::2], D[1::2]])    # deformerte koordinater
        ax = akser[r, c]                                  # ruten som skal tegnes
        ax.add_collection(PolyCollection(Xdef[IX], facecolor="yellow", edgecolor="k", linewidth=0.3))
        ax.set_xlim(Xdef[:, 0].min() - 0.5, Xdef[:, 0].max() + 0.5)   # akseutstrekning i x
        ax.set_ylim(Xdef[:, 1].min() - 0.5, Xdef[:, 1].max() + 0.5)   # akseutstrekning i y
        ax.set_aspect("equal")                            # lik skala på aksene
        ax.axis("off")                                    # fjerner akser og tall
        ax.set_title(f"x = {xv/100:.1f}, mode {m}")       # tittel med x og modenummer

plt.tight_layout()                                        # fjerner unødvendig luft
plt.savefig("modes_lbeam.png", bbox_inches="tight", pad_inches=0.05, dpi=300)   # lagrer figuren
plt.show()                                                # viser figuren