#!/usr/bin/env python3
"""
rn_korelacje.py - korelacje R_n (n-ta warstwa rownoleglego sasiedztwa)
dla lodu Ic / VII, liczone z plikow CIF (GenIce, wlasny kod MC, ...).

DEFINICJA
  R_1(X,Y): X i Y to przeciwlegle krawedzie jednego szesciokata.
  d(X,Y):   najmniejsza liczba krokow R_1 od X do Y.
  n-ta warstwa X: wszystkie Y z d(X,Y) = n, kazde liczone raz
            (w pelnej sieci 6n krawedzi: 6, 12, 18, ...).
  R_n = (zgodne - przeciwne) / (zgodne + przeciwne) po wszystkich
        parach (X, Y z n-tej warstwy X); "zgodne" = ten sam zwrot.

UZYCIE
  Kazdy folder = jeden zestaw (jedna metoda, jeden rozmiar sieci):

    python rn_korelacje.py MC/2x2x2 MC/3x3x3 ... GenIce/8x8x8 --nmax 5 -j 4

  Wyniki (domyslnie w folderze ./wyniki_Rn):
    zbiorczy.csv                   - jedna linia na (folder, n)
    <folder>_czastkowe.csv         - jedna linia na (plik, n)

  Bez podawania folderow (np. po dwukliku na plik) otwiera sie okno
  wyboru plikow: kazdy wybor to jeden zestaw, Anuluj konczy wybieranie.

    python rn_korelacje.py

  Test na wygenerowanych konfiguracjach (bez wlasnych plikow):
    python rn_korelacje.py --test

Wymagania: Python 3.8+, numpy, scipy.
"""

import argparse
import csv
import glob
import os
import re
import shlex
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import scipy.sparse as sp
from scipy.spatial import cKDTree


# ---------------------------------------------------------------------
# 1. Wczytywanie CIF (minimalny parser: komorka + petla _atom_site)
# ---------------------------------------------------------------------

def _num(s):
    """'5.123(4)' -> 5.123"""
    return float(re.sub(r"\(.*\)$", "", s))


def read_cif(path):
    """Zwraca (cell[3], typy[list], frac[N,3]). Wymaga komorki prostokatnej."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = [ln.strip() for ln in fh]
    cell = {}
    for ln in lines:
        for key in ("a", "b", "c", "alpha", "beta", "gamma"):
            if ln.startswith(f"_cell_length_{key} ") or ln.startswith(f"_cell_angle_{key} "):
                cell[key] = _num(ln.split()[1])
    for ang in ("alpha", "beta", "gamma"):
        if abs(cell.get(ang, 90.0) - 90.0) > 1e-3:
            raise ValueError(f"{path}: komorka nieprostokatna ({ang}={cell[ang]}) - nieobslugiwane")

    # znajdz petle z _atom_site_fract_x
    i = 0
    while i < len(lines):
        if lines[i].lower() == "loop_":
            j = i + 1
            cols = []
            while j < len(lines) and lines[j].startswith("_"):
                cols.append(lines[j].split()[0])
                j += 1
            if "_atom_site_fract_x" in cols:
                rows = []
                while j < len(lines) and lines[j] and not lines[j].startswith(("_", "loop_", "data_", "#")):
                    rows.extend(shlex.split(lines[j]))
                    j += 1
                nc = len(cols)
                rows = [rows[k:k + nc] for k in range(0, len(rows) - nc + 1, nc)]
                ix = {c: cols.index(c) for c in cols}
                if "_atom_site_type_symbol" in ix:
                    types = [r[ix["_atom_site_type_symbol"]] for r in rows]
                elif "_atom_site_label" in ix:
                    types = [re.match(r"[A-Za-z]+", r[ix["_atom_site_label"]]).group(0) for r in rows]
                else:
                    raise ValueError(f"{path}: brak _atom_site_type_symbol i _atom_site_label")
                frac = np.array([[_num(r[ix[f"_atom_site_fract_{a}"]]) for a in "xyz"] for r in rows])
                box = np.array([cell.get("a", 1.0), cell.get("b", 1.0), cell.get("c", 1.0)])
                return box, types, frac
            i = j
        else:
            i += 1
    raise ValueError(f"{path}: nie znaleziono petli _atom_site_fract_x/y/z")


def _is_o(t):
    t = t.strip().upper()
    return t != "OS" and re.fullmatch(r"O[0-9W_']*", t) is not None


def _is_h(t):
    t = t.strip().upper()
    return (re.fullmatch(r"[HD][0-9W_']*", t) is not None) and t not in {"HE", "HF", "HG", "HO", "HS"}


# ---------------------------------------------------------------------
# 2. Graf skierowany wiazan wodorowych
# ---------------------------------------------------------------------

def hydrogen_bond_graph(box, types, frac):
    """Zwraca (pos_O [N,3] w A, directed: dict (i,j)->True gdy wiazanie i->j,
    bonds: lista par (i<j), ostrzezenia)."""
    omask = np.array([_is_o(t) for t in types])
    hmask = np.array([_is_h(t) for t in types])
    if omask.sum() == 0 or hmask.sum() == 0:
        raise ValueError(f"brak atomow O lub H (znalezione typy: {sorted(set(types))})")
    pos_o = (np.mod(frac[omask], 1.0) * box) % box
    pos_h = (np.mod(frac[hmask], 1.0) * box) % box
    tree = cKDTree(pos_o, boxsize=box)
    _, nn = tree.query(pos_h, k=2)            # [wlasny tlen, tlen-akceptor]
    directed = set()
    for own, acc in nn:
        directed.add((int(own), int(acc)))     # kierunek: donor -> akceptor
    bonds = sorted({(min(a, b), max(a, b)) for a, b in directed})
    warn = []
    if len(bonds) != len(directed):
        warn.append(f"{len(directed) - len(bonds)} par tlenow z wiazaniem w obu kierunkach")
    deg = np.zeros(len(pos_o), int)
    for a, b in bonds:
        deg[a] += 1
        deg[b] += 1
    if np.any(deg != 4):
        warn.append(f"{int(np.sum(deg != 4))} tlenow ma liczbe wiazan rozna od 4 "
                    f"(aliasing periodyczny albo niestandardowy plik)")
    outdeg = np.zeros(len(pos_o), int)
    for a, _ in directed:
        outdeg[a] += 1
    if np.any(outdeg != 2):
        warn.append(f"{int(np.sum(outdeg != 2))} tlenow lamie regule lodu (nie 2 wodory)")
    return pos_o, directed, bonds, warn


# ---------------------------------------------------------------------
# 3. Szesciokaty i relacja R_1
# ---------------------------------------------------------------------

def find_hexagons(n_o, bonds):
    """Wszystkie 6-cykle grafu nieskierowanego jako krotki wierzcholkow."""
    nbr = [[] for _ in range(n_o)]
    for a, b in bonds:
        nbr[a].append(b)
        nbr[b].append(a)
    hexes = set()
    for a in range(n_o):
        for b in nbr[a]:
            if b < a:
                continue
            # sciezki a-b-c-d-e-f-a
            for c in nbr[b]:
                if c == a:
                    continue
                for d in nbr[c]:
                    if d in (a, b):
                        continue
                    for e in nbr[d]:
                        if e in (a, b, c):
                            continue
                        for f in nbr[e]:
                            if f in (a, b, c, d) or a not in nbr[f]:
                                continue
                            cyc = (a, b, c, d, e, f)
                            # postac kanoniczna: najmniejszy wierzcholek pierwszy, mniejszy sasiad drugi
                            k = cyc.index(min(cyc))
                            r = cyc[k:] + cyc[:k]
                            if r[1] > r[-1]:
                                r = (r[0],) + tuple(reversed(r[1:]))
                            hexes.add(r)
    return sorted(hexes)


def r1_relation(bonds, directed, hexagons):
    """Macierz znakow R_1 (+1 zgodne, -1 przeciwne) i liczba konfliktow."""
    bpos = {b: i for i, b in enumerate(bonds)}
    rows, cols, vals = [], [], []
    for v in hexagons:
        for i in range(3):
            e1 = (v[i], v[i + 1])
            e4 = (v[i + 3], v[(i + 4) % 6])
            f1 = 1 if e1 in directed else -1
            f4 = 1 if e4 in directed else -1
            # przy obiegu szesciokata przeciwlegle krawedzie sa przechodzone
            # w przeciwnych kierunkach, wiec zgodne zwroty <=> f1*f4 = -1
            s = -f1 * f4
            i1 = bpos[(min(e1), max(e1))]
            i2 = bpos[(min(e4), max(e4))]
            rows += [i1, i2]
            cols += [i2, i1]
            vals += [s, s]
    nb = len(bonds)
    msum = sp.csr_matrix((vals, (rows, cols)), shape=(nb, nb), dtype=np.int64)
    mcnt = sp.csr_matrix((np.ones(len(vals), np.int64), (rows, cols)), shape=(nb, nb))
    msum.sum_duplicates()
    mcnt.sum_duplicates()
    conflicts = int((abs(msum) != mcnt).sum()) // 2
    return msum.sign(), (mcnt > 0).astype(np.int64), conflicts


# ---------------------------------------------------------------------
# 4. Korelacje po warstwach
# ---------------------------------------------------------------------

def shell_correlations(ms, au, nmax):
    """Zwraca listy (R_n, srednia liczba krawedzi w warstwie, konflikty drog)."""
    nb = ms.shape[0]
    S = sp.identity(nb, dtype=np.int64, format="csr")
    V = S.copy()
    values, sizes, conflicts = [], [], 0
    for _ in range(nmax):
        T = S @ ms                                # suma znakow po drogach
        P = abs(S) @ au                           # liczba drog
        reach = (P > 0).astype(np.int64)
        mask = reach - reach.multiply(V)          # tylko nowe wiazania
        mask.eliminate_zeros()
        Tn = T.multiply(mask).tocsr()
        conflicts += int(((P.multiply(mask) - abs(Tn)) != 0).sum())
        S = Tn.sign().tocsr()
        V = (V + mask).tocsr()
        cnt = mask.sum()
        values.append(S.sum() / cnt if cnt else float("nan"))
        sizes.append(cnt / nb)
    return values, sizes, conflicts


def analyze_file(args):
    path, nmax = args
    try:
        box, types, frac = read_cif(path)
        pos_o, directed, bonds, warn = hydrogen_bond_graph(box, types, frac)
        hexes = find_hexagons(len(pos_o), bonds)
        ms, au, c1 = r1_relation(bonds, directed, hexes)
        values, sizes, c2 = shell_correlations(ms, au, nmax)
        if c1 or c2:
            warn.append(f"niespojne znaki: R_1 {c1}, drogi {c2}")
        return dict(path=path, n_o=len(pos_o), values=values, sizes=sizes, warn=warn, error=None)
    except Exception as exc:                      # noqa: BLE001
        return dict(path=path, error=str(exc))


# ---------------------------------------------------------------------
# 5. Program glowny
# ---------------------------------------------------------------------

def folders_to_sets(folders):
    """Folder -> zestaw (nazwa, lista plikow .cif w tym folderze)."""
    return [(f, sorted(glob.glob(os.path.join(f, "*.cif")))) for f in folders]


def set_tag(name):
    """Krotka nazwa zestawu do nazw plikow: dwa ostatnie czlony sciezki."""
    parts = [p for p in re.split(r"[\\/]+", os.path.normpath(name)) if p and not p.endswith(":")]
    return re.sub(r"[^\w.-]+", "_", "_".join(parts[-2:])).strip("_") or "zestaw"


def _f(x):
    """Liczba do CSV: zwykly zapis dziesietny, pusty napis dla NaN."""
    x = float(x)
    return "" if np.isnan(x) else f"{x:.10g}"


def run(sets, nmax, outdir, jobs):
    os.makedirs(outdir, exist_ok=True)
    summary = []
    used_tags = set()
    for folder, files in sets:
        if not files:
            print(f"[{folder}] brak plikow .cif - pomijam")
            continue
        print(f"[{folder}] {len(files)} plikow ...", flush=True)
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            res = list(ex.map(analyze_file, [(f, nmax) for f in files]))
        ok = [r for r in res if r["error"] is None]
        for r in res:
            if r["error"]:
                print(f"  BLAD {os.path.basename(r['path'])}: {r['error']}")
        nwarn = sum(1 for r in ok if r["warn"])
        if nwarn:
            first = next(r for r in ok if r["warn"])
            print(f"  ostrzezenia w {nwarn} plikach, np. {os.path.basename(first['path'])}: "
                  + "; ".join(first["warn"]))
        if not ok:
            continue
        tag = set_tag(folder)
        while tag in used_tags:
            tag += "_"
        used_tags.add(tag)
        with open(os.path.join(outdir, f"{tag}_czastkowe.csv"), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["plik", "atomow_O", "n", "R_n", "krawedzi_w_warstwie"])
            for r in ok:
                for n in range(nmax):
                    w.writerow([os.path.basename(r["path"]), r["n_o"], n + 1,
                                _f(r["values"][n]), _f(r["sizes"][n])])
        n_o = sorted({r["n_o"] for r in ok})
        vals = np.array([r["values"] for r in ok])
        sizes = np.array([r["sizes"] for r in ok])
        print(f"  atomow O: {n_o}")
        print("   n     R_n            SEM        warstwa (6n)")
        for n in range(nmax):
            v = vals[:, n][~np.isnan(vals[:, n])]
            mean = v.mean() if len(v) else float("nan")
            sem = v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else float("nan")
            shell = sizes[:, n].mean()
            full = bool(np.all(np.abs(sizes[:, n] - 6 * (n + 1)) < 1e-9))
            flag = "" if full else "  <- NIEPELNA warstwa (komorka za mala)"
            print(f"  {n + 1:2d}  {mean:+.6f}  {sem:.6f}   {shell:6.2f} ({6 * (n + 1)}){flag}")
            summary.append([folder, n_o[0] if len(n_o) == 1 else str(n_o), n + 1, _f(mean),
                            _f(sem), _f(v.std(ddof=1)) if len(v) > 1 else "",
                            len(v), _f(shell), 6 * (n + 1), int(full)])
    with open(os.path.join(outdir, "zbiorczy.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["folder", "atomow_O", "n", "R_n_srednia", "SEM", "odchylenie_std",
                    "liczba_plikow", "krawedzi_w_warstwie", "oczekiwane_6n", "warstwa_pelna"])
        w.writerows(summary)
    print(f"\nZapisano wyniki w: {os.path.abspath(outdir)}")


# ---------------------------------------------------------------------
# 6a. Wybor plikow w oknie dialogowym (tkinter)
# ---------------------------------------------------------------------

def pick_sets_gui(default_nmax):
    """Okna wyboru: kolejne zestawy plikow .cif (Anuluj konczy wybor),
    potem maksymalne n. Zwraca (zestawy, nmax, folder_wynikow)
    albo None, gdy nic nie wybrano."""
    import tkinter as tk
    from tkinter import filedialog, simpledialog, messagebox

    root = tk.Tk()
    root.withdraw()
    root.attributes("-topmost", True)
    sets = []
    start_dir = os.getcwd()
    while True:
        k = len(sets) + 1
        files = filedialog.askopenfilenames(
            parent=root,
            title=f"Zestaw {k}: zaznacz pliki .cif (Ctrl+A = wszystkie). "
                  f"Anuluj = koniec wyboru",
            initialdir=start_dir,
            filetypes=[("Pliki CIF", "*.cif"), ("Wszystkie pliki", "*.*")])
        if not files:
            break
        files = sorted(os.path.normpath(f) for f in files)
        folder = os.path.dirname(os.path.commonprefix(files) + "x")
        sets.append((folder, files))
        print(f"  zestaw {k}: {len(files)} plikow z {folder}", flush=True)
        start_dir = os.path.dirname(folder) or folder   # nastepny wybor zaczyna sie pietro wyzej
    if not sets:
        root.destroy()
        return None
    nmax = simpledialog.askinteger(
        "Maksymalny stopien n", "Do jakiego n liczyc korelacje R_n?",
        initialvalue=default_nmax, minvalue=1, maxvalue=50, parent=root)
    if nmax is None:
        nmax = default_nmax
    parents = [os.path.dirname(f) for f, _ in sets]
    base = os.path.commonpath(parents) if len({os.path.splitdrive(p)[0] for p in parents}) == 1 else parents[0]
    outdir = os.path.join(base, "wyniki_Rn")
    messagebox.showinfo(
        "Start obliczen",
        f"Wybrano zestawow: {len(sets)}\nn max = {nmax}\n\n"
        f"Wyniki zostana zapisane w:\n{outdir}\n\nPostep widac w oknie konsoli.",
        parent=root)
    root.destroy()
    return sets, nmax, outdir


# ---------------------------------------------------------------------
# 6. Test na wygenerowanych konfiguracjach lodu Ic
# ---------------------------------------------------------------------

def write_test_cif(path, L, rng, a=6.35):
    """Losowa konfiguracja lodu Ic (L x L x L komorek) spelniajaca reguly lodu."""
    M = 4 * L
    base = [(0, 0, 0), (0, 2, 2), (2, 0, 2), (2, 2, 0)]
    pts = [tuple((np.array(f) + s + 4 * np.array(c)) % M)
           for c in np.ndindex(L, L, L) for f in base for s in (0, 1)]
    idx = {p: i for i, p in enumerate(pts)}
    dirs = [np.array(d) for d in np.ndindex(2, 2, 2)]
    dirs = [2 * d - 1 for d in dirs]
    nbr = {p: [tuple((np.array(p) + d) % M) for d in dirs if tuple((np.array(p) + d) % M) in idx]
           for p in pts}
    out = {p: [q for q in nbr[p] if ((np.array(q) - p + M // 2) % M - M // 2)[2] > 0] for p in pts}
    for _ in range(20 * len(pts)):                # losowe odwracanie petli
        cur = pts[rng.integers(len(pts))]
        walk, seen = [cur], {cur: 0}
        while True:
            nxt = out[cur][rng.integers(2)]
            if nxt in seen:
                cyc = walk[seen[nxt]:] + [nxt]
                break
            seen[nxt] = len(walk)
            walk.append(nxt)
            cur = nxt
        for u, v in zip(cyc, cyc[1:]):
            out[u].remove(v)
            out[v].append(u)
    with open(path, "w") as fh:
        size = a * L
        fh.write(f"data_test\n_cell_length_a {size}\n_cell_length_b {size}\n_cell_length_c {size}\n"
                 "_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n"
                 "loop_\n_atom_site_label\n_atom_site_type_symbol\n"
                 "_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n")
        k = 0
        for p in pts:
            fh.write(f"O{idx[p]} O {p[0] / M:.6f} {p[1] / M:.6f} {p[2] / M:.6f}\n")
        for p in pts:
            for q in out[p]:
                d = (np.array(q) - p + M // 2) % M - M // 2
                h = (np.array(p) + 0.35 * d) / M % 1.0
                fh.write(f"H{k} H {h[0]:.6f} {h[1]:.6f} {h[2]:.6f}\n")
                k += 1


def self_test():
    import tempfile
    rng = np.random.default_rng(0)
    tmp = tempfile.mkdtemp(prefix="rn_test_")
    folders = []
    for L in (2, 3, 4, 6):
        d = os.path.join(tmp, f"{L}x{L}x{L}")
        os.makedirs(d)
        for k in range(4):
            write_test_cif(os.path.join(d, f"test_{k}.cif"), L, rng)
        folders.append(d)
    run(folders_to_sets(folders), 5, os.path.join(tmp, "wyniki"), jobs=1)
    print("\nOczekiwane: warstwy 6/12/18/24/30 pelne dla n <= L-1, brak ostrzezen o znakach.")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("foldery", nargs="*", help="foldery z plikami .cif (jeden folder = jeden zestaw); "
                                               "bez folderow otwiera sie okno wyboru plikow")
    ap.add_argument("--nmax", type=int, default=5, help="maksymalny stopien n (domyslnie 5)")
    ap.add_argument("-o", "--out", default=None, help="folder na wyniki (domyslnie wyniki_Rn)")
    ap.add_argument("-j", "--jobs", type=int, default=os.cpu_count(), help="liczba procesow")
    ap.add_argument("--test", action="store_true", help="test na wygenerowanych konfiguracjach")
    a = ap.parse_args()
    if a.test:
        self_test()
    elif a.foldery:
        run(folders_to_sets(a.foldery), a.nmax, a.out or "wyniki_Rn", a.jobs)
    else:
        print("Wybierz pliki .cif w oknie dialogowym: jeden wybor = jeden zestaw "
              "(jedna metoda, jeden rozmiar). Anuluj konczy wybieranie.", flush=True)
        try:
            picked = pick_sets_gui(a.nmax)
        except Exception as exc:                  # noqa: BLE001
            print(f"Nie udalo sie otworzyc okna wyboru ({exc}). "
                  f"Podaj foldery w wierszu polecen, np.: python rn_korelacje.py MC/4x4x4 GenIce/4x4x4")
            picked = None
        if picked:
            sets, nmax, outdir = picked
            run(sets, nmax, a.out or outdir, a.jobs)
        else:
            print("Nie wybrano zadnych plikow.")
        try:
            input("\nNacisnij Enter, aby zamknac...")
        except EOFError:
            pass


if __name__ == "__main__":
    main()
