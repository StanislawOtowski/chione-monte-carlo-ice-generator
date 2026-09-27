#!/usr/bin/env python3
"""
porownanie_metod.py - statystyczne porownanie korelacji C_n z dwoch metod
(Monte Carlo i GenIce) na podstawie plikow zbiorczy.csv z rn_korelacje.py.

Dla kazdego rozmiaru sieci i stopnia n (tylko pelne warstwy):
    z_n = (C_MC - C_GI) / sqrt(u_MC^2 + u_GI^2)            (wzor eq:z)
a dla kazdego rozmiaru laczny test:
    chi^2 = sum_n z_n^2,  stopnie swobody = liczba skladnikow,
    p = P(chi^2_k >= chi^2)                                  (wzor eq:chi2)
Uwaga: C_n roznych stopni liczone z tych samych konfiguracji sa skorelowane,
wiec p z testu lacznego jest przyblizone.

UZYCIE
    python porownanie_metod.py zbiorczy.csv [poprawki.csv ...] [--tex tabela.tex]

Mozna podac kilka plikow CSV: wiersze z pozniejszych plikow zastepuja wiersze
o tym samym folderze i n z plikow wczesniejszych (np. przeliczony folder).
Metoda jest rozpoznawana po nazwie folderu (domyslnie: zawiera "MonteCarlo"
albo "GenIce"; zmiana opcjami --mc i --gi).

Wymagania: Python 3.8+, scipy.
"""

import argparse
import csv
import math
import re
import sys

from scipy.stats import chi2, norm


def load(paths):
    """Zwraca slownik (folder, n) -> wiersz; pozniejsze pliki nadpisuja wczesniejsze."""
    rows = {}
    for path in paths:
        with open(path, newline="", encoding="utf-8-sig") as fh:
            for r in csv.DictReader(fh):
                rows[(r["folder"], int(r["n"]))] = r
    return rows


def collect(rows, mc_pat, gi_pat):
    """data[metoda][liczba_atomow][n] = (srednia, SEM, liczba_plikow)."""
    data = {"MC": {}, "GI": {}}
    for (folder, n), r in rows.items():
        if r.get("warstwa_pelna", "1").strip() != "1":
            continue
        if re.search(mc_pat, folder, re.I):
            m = "MC"
        elif re.search(gi_pat, folder, re.I):
            m = "GI"
        else:
            print(f"UWAGA: nie rozpoznano metody dla folderu {folder} - pomijam", file=sys.stderr)
            continue
        try:
            atoms = int(r["atomow_O"])
        except ValueError:
            print(f"UWAGA: folder {folder} zawiera pliki o roznych rozmiarach "
                  f"({r['atomow_O']}) - pomijam", file=sys.stderr)
            continue
        if not r["SEM"]:
            continue
        data[m].setdefault(atoms, {})[n] = (float(r["R_n_srednia"]), float(r["SEM"]),
                                           int(r["liczba_plikow"]))
    return data


def compare(data):
    """Lista wynikow dla rozmiarow wystepujacych w obu metodach."""
    out = []
    for atoms in sorted(set(data["MC"]) & set(data["GI"])):
        ns = sorted(set(data["MC"][atoms]) & set(data["GI"][atoms]))
        zs = []
        for n in ns:
            a, ua, _ = data["MC"][atoms][n]
            b, ub, _ = data["GI"][atoms][n]
            zs.append((n, (a - b) / math.hypot(ua, ub), a - b))
        x2 = sum(z * z for _, z, _ in zs)
        p = chi2.sf(x2, len(zs)) if zs else float("nan")
        L = round(atoms ** (1 / 3) / 2)          # 8 atomow na komorke elementarna Ic
        out.append(dict(atoms=atoms, L=L, zs=zs, chi2=x2, dof=len(zs), p=p))
    return out


def p_tex(p):
    """p w postaci m\\cdot10^{e} z jedna cyfra znaczaca (zaokraglenie w gore od polowy)."""
    e = math.floor(math.log10(p))
    m = math.floor(p / 10 ** e + 0.5)
    if m == 10:
        m, e = 1, e + 1
    return f"${m}\\cdot10^{{{e}}}$"


def z_tex(z):
    return f"${z:.1f}$".replace(".", "{,}")


def write_tex(res, path, nmax=5):
    lines = []
    for i, r in enumerate(res):
        zmap = {n: z for n, z, _ in r["zs"]}
        cells = [z_tex(zmap[n]) if n in zmap else "" for n in range(1, nmax + 1)]
        end = r" \\" if i < len(res) - 1 else ""
        lines.append(f"        ${r['L']}\\times{r['L']}\\times{r['L']}$ & {r['atoms']} & "
                     + " & ".join(cells) + f" & {p_tex(r['p'])}" + end)
    zcols = " & ".join(f"$z_{n}$" for n in range(1, nmax + 1))
    tex = r"""% Tabela porownania metod MC i GenIce -- wczytywana przez \input{porownanie_z_tabela}
% z_n = (C_n^MC - C_n^GI) / sqrt(u_MC^2 + u_GI^2)  (wzor eq:z)
% p   = prawdopodobienstwo z lacznego testu chi^2 = sum_n z_n^2 (wzor eq:chi2),
%       liczba stopni swobody = liczba pelnych warstw dla danego rozmiaru.
% Wartosci policzone z tych samych danych co mc_tabela.tex i genice_tabela.tex.
\begin{table}[htbp]
    \centering
    \small
    \caption{\label{tab:porownanie_z}Różnice korelacji $\bar{C}_n$ otrzymanych metodą Monte Carlo i programem GenIce 2, wyrażone wzorem~\eqref{eq:z}, oraz prawdopodobieństwo $p$ przypadkowej zgodności obu metod według łącznego testu~\eqref{eq:chi2}. Wartości dodatnie oznaczają, że korelacja z GenIce 2 jest niższa (bardziej ujemna).}
    \begin{tabular}{c|c|""" + "c" * nmax + r"""|c}
        rozmiar & liczba atomów tlenu & """ + zcols + r""" & $p$ \\
        \hline\hline
""" + "\n".join(lines) + r"""
    \end{tabular}
\end{table}
"""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(tex)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", nargs="+", help="pliki zbiorczy.csv (pozniejsze nadpisuja wczesniejsze)")
    ap.add_argument("--mc", default=r"monte\s*carlo", help="wzorzec nazwy folderu MC (regex)")
    ap.add_argument("--gi", default=r"genice", help="wzorzec nazwy folderu GenIce (regex)")
    ap.add_argument("--tex", help="zapisz tabele LaTeX do tego pliku")
    ap.add_argument("--csv-out", help="zapisz wyniki (rozmiar, n, roznica, z) do pliku CSV")
    a = ap.parse_args()

    data = collect(load(a.csv), a.mc, a.gi)
    res = compare(data)
    if not res:
        sys.exit("Brak rozmiarow wspolnych dla obu metod.")

    print("rozmiar  atomy   z_n (n=1..)                              chi2   k   p")
    for r in res:
        zs = "  ".join(f"{z:+5.1f}" for _, z, _ in r["zs"])
        print(f"{r['L']}x{r['L']}x{r['L']}  {r['atoms']:5d}   {zs:40s} {r['chi2']:6.1f}  {r['dof']}  {r['p']:.3g}")

    allz = [z for r in res for _, z, _ in r["zs"]]
    pos = sum(z > 0 for z in allz)
    print(f"\nliczba porownan: {len(allz)};  z > 0: {pos};  |z| >= 2: {sum(abs(z) >= 2 for z in allz)};"
          f"  |z| >= 3: {sum(abs(z) >= 3 for z in allz)}")
    print(f"dla porownania: P(|z| >= 2) = {2 * norm.sf(2):.3f},  P(|z| >= 3) = {2 * norm.sf(3):.4f}"
          f" przy zgodnosci metod")

    if a.tex:
        write_tex(res, a.tex)
        print(f"Zapisano tabele LaTeX: {a.tex}")
    if a.csv_out:
        with open(a.csv_out, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["atomow_O", "L", "n", "roznica_MC_minus_GI", "z"])
            for r in res:
                for n, z, d in r["zs"]:
                    w.writerow([r["atoms"], r["L"], n, f"{d:.10g}", f"{z:.6f}"])
            w.writerow([])
            w.writerow(["atomow_O", "L", "chi2", "stopnie_swobody", "p"])
            for r in res:
                w.writerow([r["atoms"], r["L"], f"{r['chi2']:.6f}", r["dof"], f"{r['p']:.6g}"])
        print(f"Zapisano wyniki CSV: {a.csv_out}")


if __name__ == "__main__":
    main()
