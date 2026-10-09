# barometer

Analyse exhaustive de biomarqueurs (BMK) sur les données d'édition RAIN.

Le pipeline lit les tables TSV produites par RAIN (aggregates / features / sites),
découpe chaque table en **sections** (par type de feature, type de compte, mode
d'isoforme…), et pour chaque section exécute une batterie d'analyses
statistiques et de machine learning. Les résultats sont agrégés en un **classement
global** par type de valeur, puis un **second passage** est lancé sur les BMK
significatifs. Un rapport HTML interactif est généré à la fin.

## Structure du package

| Fichier | Rôle |
|---------|------|
| `__main__.py` | Point d'entrée CLI (`python -m barometer`), orchestration, parallélisation, retry, `--resume`, second passage, merge. |
| `analysis.py` | Les analyses par section (QC, batch, descriptive, multivarié, différentiel, corrélation, ranking, classification, stabilité, heatmap). |
| `ranking.py` | Classement global par mtype + diagrammes de significativité. |
| `merge.py` | Fusion des classements globaux par mtype en un classement « vraiment global » + second passage. |
| `report.py` | Génération du rapport HTML interactif. |
| `utils.py` | Chargement de données, création de UID, parsing des colonnes d'échantillons, helpers. |
| `barometer_beta_binomial.R` | Script R (glmmTMB) pour le mode beta-binomial. |
| `tests/` | Tests synthétiques (voir [Tests](#tests)). |

## Déroulé du script

### 1. Chargement des données

Les fichiers TSV (`--aggregates`, `--features`, `--sites`) sont lus, leurs
colonnes harmonisées, un **UID** est créé pour identifier chaque ligne, et les
colonnes à faible cardinalité (`Mtype`, `Ptype`, `Type`, `Ctype`, `Mode`,
`Strand`, `SeqID`) sont converties en `category` pour économiser de la mémoire.

Les colonnes d'échantillons sont détectées automatiquement
(`parse_sample_columns`). Chaque **type de valeur** (`espf`, `espr`, …) est
analysé séparément.

### 2. Découpage en sections

Pour chaque type de valeur, les lignes sont groupées en sections selon les
attributs de feature (Ptype × Ctype × Mode, etc.). Chaque section devient une
**tâche** indépendante, ce qui permet la parallélisation.

### 3. Analyse par section (`analyze_section`)

Chaque section exécute, dans l'ordre, les étapes suivantes (chacune dans son
sous-dossier numéroté) :

| # | Étape | Dossier | Portée |
|---|-------|---------|--------|
| 0 | `filter1` — suppression des BMK à variance quasi nulle | — | tous les BMK |
| 0 | `filter2` — limitation à 10 000 BMK (top variance/trials) | `data.csv` | tous les BMK |
| 1 | `qc` — contrôle qualité | `1_qc/` | tous les BMK |
| 2 | `batch` — détection d'effet de lot | `2_batch/` | tous les BMK |
| 3 | `descriptive` — statistiques descriptives | `3_descriptive/` | tous les BMK |
| 4 | `differential` — tests différentiels | `5_differential/` | tous les BMK |
| — | **Sélection** des BMK différentiels (`_select_differential_bmks`) | — | — |
| 5 | `multivariate` — PCA, clustering | `4_multivariate/` | BMK sélectionnés |
| 6 | `correlation` — réseau de corrélation | `6_correlation/` | BMK sélectionnés |
| 7 | `ranking` — sélection de features / classement | `7_ranking/` | BMK sélectionnés |
| 8 | `classification` — modélisation prédictive | `8_classification/` | BMK sélectionnés |
| 9 | `stability` — concordance des réplicats | `9_stability/` | BMK sélectionnés |
| 10 | `heatmap` — carte de chaleur | (racine) | BMK sélectionnés |

**Conception clé** : les étapes `qc`, `batch`, `descriptive` et `differential`
sont faites sur **tous** les BMK de la section. Ensuite, les BMK différentiels
sont sélectionnés, et les étapes qui en dépendent (`multivariate`, `correlation`,
`ranking`, `classification`, `stability`, `heatmap`) ne tournent que sur ces BMK
sélectionnés. Cela évite de gaspiller du calcul sur des BMK non informatifs.

Le paramètre `enabled_tests` permet de restreindre les étapes (utilisé par le
second passage, qui saute `qc`/`batch`).

### 4. Tests différentiels

Le pipeline calcule **toujours tous** les tests pour chaque BMK :

- **Tests globaux** (tous les groupes) : Kruskal-Wallis, ANOVA, Welch ANOVA.
- **Tests par paires** (X vs Y) : Mann-Whitney U, Student t-test, Welch t-test.
- **Tests d'hypothèses** (mode `auto`) : Shapiro-Wilk, Bartlett.

Tous reçoivent une correction FDR (Benjamini-Hochberg) → colonnes `*_padj`.
Résultats dans `5_differential/differential_results.csv`.

**Mode beta-binomial** (`--stat-test beta-binomial`) : au lieu de tester les
proportions arrondies, un GLMM beta-binomial (R/glmmTMB) est ajusté sur les
**comptes bruts** DRIP (`{sample}::successes` / `{sample}::trials`). Le paramètre
de dispersion absorbe la surdispersion extra-binomiale typique du séquençage.
Deux tests de Wald en découlent : un **global** (1/BMK) et des **par paires**
(C(n,2)/BMK). Nécessite `Rscript` + `glmmTMB` (conteneur Barometer).

### 5. Parallélisation, retry et `--resume`

Les sections sont exécutées en parallèle via `ProcessPoolExecutor`
(`-j N`, ou `-1` pour tous les cœurs). Un **superviseur** gère :

- **Timeout par tâche** : chaque tâche a un budget `SIGALRM` (300 s par défaut).
- **Retry avec timeout doublé** : si une tâche échoue pour une raison liée au
  temps (timeout, slow, deadlock, orphaned), elle est retentée jusqu'à 3 fois
  avec un timeout qui double : 300 s → 600 s → 1200 s → 2400 s. Les erreurs de
  code (crash, exception) ne sont **pas** retentées.
- **Watchdog** : détection des tâches lentes et des deadlocks (aucune
  progression).
- **`failed_analyses.tsv`** : à la fin, toutes les sections qui ont échoué sont
  listées dans `<outdir>/failed_analyses.tsv` (colonnes : `vtype`, `mtype`,
  `section_key`, `section_name`, `section_outdir`, `reason`, `detail`,
  `timestamp`).

**`--resume`** : pour relancer uniquement les sections en échec :

```bash
python -m barometer ... --resume <outdir>/failed_analyses.tsv
```

Seules les sections du TSV sont relancées. Si certaines échouent à nouveau, elles
figurent dans le **nouveau** `failed_analyses.tsv` → on peut itérer `--resume`
jusqu'à ce que le fichier soit vide.

### 6. Classement global et second passage

Pour chaque type de valeur et chaque mtype (aggregate / feature / sites),
`global_ranking` agrège les tables différentielles des sections en un classement
unique, et identifie les BMK significatifs :

- `significant_bmks_any_comparison/` — union de toutes les comparaisons.
- `significant_bmks_global_comparison/` — test global (omnibus) seulement.
- `significant_bmks_pairwise_comparison/` — par comparaison de paires.

Un **second passage** (`analyze_section` avec `enabled_tests` restreint, sans
`qc`/`batch`) est alors lancé sur ces BMK significatifs, pour les ré-analyser
dans leur contexte propre.

### 7. Merge (classement vraiment global)

`--merge` combine les classements globaux par mtype (de plusieurs runs) en un
classement **vraiment global**, puis relance un second passage par redistribution
par type de valeur.

### 8. Rapport HTML

`report.py` génère un rapport HTML interactif à partir des dossiers de résultats
(`-r/--results`, ou `--espf`/`--espr` par mtype, et `--merged` pour le merge).

## Utilisation

```bash
# Analyse par défaut (cascade padj + pval, 500 BMK max) :
python -m barometer -a data.tsv -o results/ -j 4

# Seulement Kruskal-Wallis :
python -m barometer -a data.tsv -o results/ --bmk-filter kruskal_padj

# Welch en priorité, fallback Kruskal, 1000 BMK :
python -m barometer -a data.tsv -o results/ \
    --stat-test welch --bmk-filter welch_padj kruskal_padj --max-bmks 1000

# Beta-binomial sur les comptes bruts DRIP (conteneur Barometer) :
python -m barometer -a data.tsv -o results/ --stat-test beta-binomial

# Relancer les sections en échec :
python -m barometer -a data.tsv -o results/ --resume results/failed_analyses.tsv

# Fusionner en classement vraiment global :
python -m barometer --merge --results-dir run1/ run2/ -o merged/ \
    --raw-inputs '{"espf": {"aggregates": "a.tsv", "features": "f.tsv"}}'

# Générer le rapport HTML :
python -m barometer.report -r results/ -o report.html
```

### Options principales

| Option | Description |
|--------|-------------|
| `-a/--aggregates`, `-f/--features`, `--sites` | Fichiers TSV d'entrée (au moins un). |
| `-o/--outdir` | Dossier de sortie (défaut : `barometer_results`). |
| `-j/--jobs` | Nombre de jobs parallèles (`-1` = tous les cœurs). |
| `-v/--value-types` | Types de valeur à analyser (ex. `espf espr`). |
| `--agg-levels` | Niveaux d'agrégation : `global`, `sequence`, `feature`. |
| `--feature-types` | Types de feature à analyser. |
| `--stat-test` | `auto`, `parametric`, `nonparametric`, `welch`, `kruskal`, `beta-binomial`. |
| `--bmk-filter` | Liste prioritaire de colonnes `*_padj` pour la sélection des BMK. |
| `--max-bmks` | Nombre max de BMK pour les modèles ML (défaut : 500). |
| `--merge` | Mode fusion des classements globaux. |
| `--results-dir` | Dossiers de résultats pour `--merge`. |
| `--raw-inputs` | JSON vtype → TSV bruts (pour le second passage du merge). |
| `--resume` | Fichier `failed_analyses.tsv` à relancer. |

## Tests

Le dossier `tests/` contient des **tests synthétiques** qui valident la logique
du superviseur sans dépendre des bibliothèques lourdes (seaborn, sklearn,
statsmodels) ni de données réelles. Ils réimplémentent verbatim les morceaux de
code concernés et les pilotent avec des workers factices.

| Test | Ce qu'il vérifie |
|------|------------------|
| `test_supervisor_retry.py` | Le retry avec timeout doublé : une tâche qui échoue 1× puis réussit est retentée une fois ; une tâche qui échoue toujours est retentée 3× (600/1200/2400 s) puis consignée dans `failed_analyses` ; une tâche qui réussit immédiatement n'est pas retentée. |
| `test_resume_filter.py` | Le mode `--resume` : le chargement du TSV et le filtrage des tâches ne conservent que les sections listées comme en échec. |

Ces tests n'importent pas le package `barometer` (qui tirerait seaborn/sklearn/
statsmodels) : ils sont donc exécutables avec un `python3` minimal (pandas +
multiprocessing).

### Cible Python 3.11+ et gates de version

Le code de production (`__main__.py`) cible **Python 3.11+** (le container
barometer tourne sous **3.12**). Deux points en dépendent :

- **`max_tasks_per_child`** (3.11+) : le superviseur recycle les workers
  (`ProcessPoolExecutor(max_workers=n_jobs, max_tasks_per_child=20, ...)`).
- **Unification de `TimeoutError`** (3.11+) : `concurrent.futures.TimeoutError`
  **est** le `TimeoutError` builtin, donc un simple `except TimeoutError:`
  suffit (pas de double-catch).

`test_supervisor_retry.py` reste **exécutable en 3.9 local** (env de dev) grâce
à des gates de version :

```python
_TIMEOUT_EXC = (TimeoutError,) if sys.version_info >= (3, 11) else (TimeoutError, FuturesTimeoutError)
_pool_kwargs = {"max_tasks_per_child": 20} if sys.version_info >= (3, 11) else {}
```

En 3.11+ le test reproduit fidèlement le vrai superviseur (`max_tasks_per_child=20`,
`except TimeoutError:`) ; en 3.9 il dégrade proprement ces deux points.

```bash
cd bin/barometer
python3 tests/test_supervisor_retry.py
python3 tests/test_resume_filter.py
```

Chaque test affiche `TEST PASSED: ...` en cas de succès et lève une
`AssertionError` sinon.
