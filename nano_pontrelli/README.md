# Nano-contenitori di Pontrelli — Progetto di Modelli Fisiologici

Simulazione del rilascio di farmaci da nano-contenitori sensibili al pH, basata sul modello di Pontrelli.

## Contenuto del repository

| File | Descrizione |
|------|-------------|
| `progetto3d_0401.mph` | Simulazione COMSOL 3D (versione del 04/01) |
| `progettoR_2901.mph` | Simulazione COMSOL assialsimmetrica (versione del 29/01) |
| `cleanind_data.ipynb` | Notebook Python per la ripulitura e l'analisi dei dati |
| `massa_rilasciata2901.m` | Script MATLAB con i grafici della massa rilasciata |

## Descrizione

Il progetto studia la diffusione e il rilascio controllato di una molecola farmacologica da nano-contenitori polimerici sensibili al pH. Le simulazioni COMSOL risolvono le equazioni di diffusione-reazione nel dominio della nanoparticella, mentre il notebook Python elabora i dati esportati e lo script MATLAB produce i grafici finali.

## Dipendenze

- **COMSOL Multiphysics** >= 6.x (per aprire i file `.mph`)
- **Python** >= 3.9 con `numpy`, `pandas`, `matplotlib` (per il notebook)
- **MATLAB** >= R2022a (per lo script `.m`)
