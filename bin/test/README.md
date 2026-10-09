# Tests

Tests unitaires des scripts d'analyse de `bin/` (drip, barometer, pluviometer).

## Contenu

| Fichier | Cadre | Ce qu'il teste | Environnement |
|---------|-------|----------------|---------------|
| `test_drip_counts.py` | unittest | `drip.merge_samples` : ratios espf/espr, colonnes successes/trials, seuil de couverture, `preserve_covered_features` (ajoute les colonnes de compte sans modifier le filtrage de lignes), filtres de couverture `min_samples` (AND avec `min_samples_pct` sur le compteur global) et `min_group_samples` (AND avec `min_group_samples_pct` sur le même groupe, protège du bruit 1-sample), filtres d'édition `min_group_samples_edited` / `min_group_samples_pct_edited` (AND sur le même groupe, défaut 1/0), zéros couverts comptent pour la couverture mais pas pour l'édition | Conteneur `drip` (Python 3.12, numpy, pandas) |
| `test_barometer_beta_binomial.py` | unittest | Pont `barometer_analyze.run_beta_binomial_analysis` vers le script R `barometer_beta_binomial.R` (mapping des résultats, agrégation des répliques techniques) | Conteneur `barometer` (Python 3.12, seaborn, R 4.4 + glmmTMB) |
| `test_site_format.py` | pytest | Format de sortie 16 colonnes de `pluviometer.site_file_writer` et traitement en aval par `drip.py` (lancé en subprocess) | Conteneur `pluviometer` (Python 3.12, Biopython) ; le subprocess `drip.py` n'a besoin que de pandas/numpy (conteneur `drip`) |

## Exécution

Depuis la racine du dépôt :

```bash
# Tous les tests (pytest collecte aussi les classes unittest)
python -m pytest bin/test/ -v

# Ou avec unittest uniquement
python -m unittest discover -s bin/test -v
```

Un seul fichier :

```bash
python -m pytest bin/test/test_drip_counts.py -v
```

## Environnements

Les tests sont conçus pour s'exécuter dans les conteneurs du pipeline (voir `containers/`) :

| Conteneur | Dépendances clés | Construit depuis |
|-----------|------------------|------------------|
| `drip` | Python 3.12, numpy, pandas | `containers/docker/drip/Dockerfile` |
| `barometer` | Python 3.12, seaborn, R 4.4 + glmmTMB | `containers/docker/barometer/Dockerfile` |
| `pluviometer` | Python 3.12, Biopython, pysam | `containers/docker/pluviometer/Dockerfile` |

Exemple d'exécution dans un conteneur :

```bash
docker run --rm -v "$PWD":/app -w /app/bin <image-drip> \
    python -m pytest test/ -v
```

## Notes

- Les tests ajoutent `bin/` au `sys.path` eux-mêmes (via `Path(__file__).resolve().parent.parent`) : aucun `PYTHONPATH` ni `conftest.py` n'est nécessaire.
- `test_barometer_beta_binomial.py` contient un test optionnel (`test_glmmtmb_detects_synthetic_condition_effect`) qui est **sauté automatiquement** si `Rscript` n'est pas disponible (disponible dans le conteneur Barometer).
- `test_site_format.py` stubbe `Bio.SeqFeature` pour pouvoir s'exécuter sans BioPython installé.
- Les tests du package pluviometer (`pluviometer/test_site_counts.py`) restent dans le package et s'exécutent séparément :
  ```bash
  cd bin && python -m unittest pluviometer.test_site_counts -v
  ```


### quickstats.sh
Quick statistics tool for variant calling output files. Counts A→G and C→T edits per contig/chromosome from Reditools or Jacusa2 output files.

**Usage:**
```bash
./quickstats.sh -i INPUT_FILE -f FORMAT
# FORMAT: "jacusa2" or "reditools"
``` 