# Monte Carlo konfiguracji protonowych lodu Ic i VII

[English version](README.md)

Symulacja Monte Carlo (algorytm Metropolisa) rozmieszczenia wodorów w sieci lodu Ic i lodu VII ze stopniowym schładzaniem. Energia to kara za łamanie reguł lodu: przy każdym tlenie powinny leżeć dokładnie 2 wodory. Wynikiem są pliki CIF z końcową konfiguracją.

## Instalacja

```
pip install -r requirements.txt
```

## Uruchomienie

```
python main.py
```

Parametry ustawia się w `main.py`:

- w bloku `if __name__=='__main__':` liczbę rdzeni (`Pool(4)`) i listę rozmiarów sieci (`[4]`), które
  liczą się równolegle, po jednym zadaniu na każdy rozmiar;
- w argumentach `generate_and_iterate_net`:
  - `ice_type`: `'1c'` albo `'7'`,
  - `temperatures`: lista kT (od najwyższej do najniższej),
  - `iterations_at_one_temperature_at_stable_state`: liczba próbek na temperaturę
    (przed nimi wykonuje się połowę tej liczby sweepów na dojście do równowagi),
  - `cif_file_generation`: czy zapisywać CIF.

Rozmiar sieci `net_size = n` oznacza n×n×n komórek elementarnych: 8n³ tlenów dla lodu Ic i 16n³ dla lodu VII (dwie przenikające się podsieci Ic).

## Wyniki

- Na konsolę, dla każdej temperatury: energia, momenty energii na cząsteczkę (⟨e⟩…⟨e⁴⟩) i ciepło właściwe cV na cząsteczkę.
- Pliki CIF trafiają do `Plotting/CIF/`, nazwa `Ice_<typ>_<liczba tlenów>_oxygens_<data>_<godzina>.<losowe>.cif`.
- Na koniec otwiera się wykres 3D struktury (`plt.show()` blokuje, dopóki okno nie zostanie zamknięte).

## Struktura

| plik | rola |
|---|---|
| `main.py` | uruchomienie: generacja sieci, MC, zapis CIF, wykres |
| `Generating/Ic_ice_net_generator.py` | generator sieci Ic z losową konfiguracją wodorów |
| `Ice_Structure_Changes/change_hydrogen_configuration.py` | pętla MC (Metropolis) po temperaturach |
| `Analyzys/energy_calculations.py` | energia całkowita i lokalna zmiana energii |
| `Plotting/saving_to_CIF_file.py` | zapis CIF na podstawie szablonów `Ice_1c_template.txt`, `Ice_7_template.txt` |
| `Plotting/Ic_and_7_ice_net_plotter.py` | wykres 3D |
