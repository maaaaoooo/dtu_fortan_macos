import os                                      # filstier
import subprocess                              # kjører Fortran-programmet
import numpy as np                             # arrays og linspace
import matplotlib.pyplot as plt                # plotting

h = 0.25                                       # elementstørrelse (prøv 0.125 for konvergenssjekk)
n = round(1 / h)                               # antall elementer per lengdeenhet
E, nu, rho = 50.0, 0.0, 2.0                    # materialdata fra oppgaven
FIL = "lbeam"                                  # navn på inputfilen
MAPPE = "."                                    # mappen der main_program ligger
EPS = 1e-3                                     # minste tykkelse, unngår t = 0 (singulær K)

def lag_input(x):                              # skriver inputfilen for én verdi av x
    t_soyle = max(x, EPS)                      # søyle: t = x
    t_bjelke = max(1 - x, EPS)                 # bjelke: t = 1 - x
    soyle = [(i, j) for i in range(1 * n) for j in range(5 * n)]          # celler i søylen (x 0-1, y 0-5)
    bjelke = [(i, j) for i in range(10 * n) for j in range(5 * n, 6 * n)] # celler i bjelken (x 0-10, y 5-6)
    noder = {}                                 # (i, j) i rutenettet -> nodenummer

    def nodenr(i, j):                          # gir nodenummer, lager nytt første gang
        if (i, j) not in noder:                # hvis punktet ikke har nummer ennå
            noder[(i, j)] = len(noder) + 1     # neste ledige nummer, uten hull
        return noder[(i, j)]                   # returnerer nummeret

    def celle(i, j):                           # fire hjørner mot klokka
        return [nodenr(i, j), nodenr(i + 1, j), nodenr(i + 1, j + 1), nodenr(i, j + 1)]

    el_soyle = [celle(i, j) for i, j in soyle]     # elementer i søylen
    el_bjelke = [celle(i, j) for i, j in bjelke]   # elementer i bjelken (deler noder i hjørnet)

    L = ["/TITLE, L-FRAME", "/PREP7", "ET,1,PLANE42"]            # tittel, preprosessor, elementtype
    for m, t in [(1, t_soyle), (2, t_bjelke)]:                   # materiale 1 = søyle, 2 = bjelke
        L += [f"MP,EX,{m}, {E}", f"MP,PRXY,{m}, {nu}", f"MP,DENS,{m}, {rho}", f"R,{m}, {t:.8f}"]
    for (i, j), nr in sorted(noder.items(), key=lambda kv: kv[1]):  # noder i nummerrekkefølge
        L.append(f"N,{nr}, {i * h:.8f}, {j * h:.8f}, 0.")           # nodekoordinater
    enr = 0                                    # elementteller, må gå 1, 2, 3, ... uten hull
    for m, elementer in [(1, el_soyle), (2, el_bjelke)]:         # gruppér elementer per materiale
        L.append(f"MAT,{m}")                   # velger materiale for EN-linjene under
        for nodene in elementer:               # ett element om gangen
            enr += 1                           # neste elementnummer
            L.append(f"EN,{enr}, " + ", ".join(map(str, nodene)))   # elementets nodeliste
    fast = sorted(nr for (i, j), nr in noder.items() if j == 0 or i == 10 * n)  # bunn av søyle og høyre side av bjelke
    for nr in fast:                            # fastholdte noder
        L += [f"D,{nr}, UX, 0.0", f"D,{nr}, UY, 0.0"]            # begge retninger låst
    L += ["FINISH", "/SOLU", "antype, modal", "FINISH"]          # modalanalyse
    with open(os.path.join(MAPPE, FIL), "w") as f:               # skriver filen
        f.write("\n".join(L) + "\n")

def les_omega(ut, antall=6):                   # henter omega ut av programmets utskrift
    omega = []                                 # liste for omega-verdiene
    for linje in ut.splitlines():              # går gjennom hver linje
        d = linje.split()                      # deler linjen i ord
        if len(d) == 4 and d[0].isdigit():     # linjene fra eigen: modenr, lambda, omega, f
            omega.append(float(d[2]))          # tredje tall er omega
    return omega[:antall]                      # de første modene

xs = np.linspace(0, 1, 20)                     # 20 verdier av x fra 0 til 1
omega = np.zeros((len(xs), 6))                 # rad = x-verdi, kolonne = modus
for k, x in enumerate(xs):                     # løkke over x
    lag_input(x)                               # lager inputfilen for denne x
    ut = subprocess.run(["./main_program"], input=FIL + "\n", capture_output=True,
                        text=True, cwd=MAPPE).stdout   # kjører programmet, fanger utskriften
    omega[k, :] = les_omega(ut)                # lagrer de seks frekvensene i rad k
np.save("omega_lbeam.npy", omega)              # lagrer arrayet til fil

for i in range(6):                             # én kurve per modus
    plt.plot(xs, omega[:, i], marker="^", alpha=1, label=f"Mode {i + 1}")
plt.xlabel("x")                                # x-aksen
plt.ylabel(r"$\omega$ [rad/s]")                # y-aksen (del på $2\pi$ for Hz)
plt.legend(fontsize=14, markerscale=0.5, handlelength=1.5, frameon=True)         # tegnforklaring
plt.grid()
plt.show()                                     # viser plottet