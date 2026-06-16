# BCR Lineage Tracer (B cell Linage Tree Maker)
In this project, I intend to create a Python tool for making B cell receptor (BCR) clonal lineage trees from single-cell sequencing data.
The tool infers the ancestral germline BCR, traces somatic hypermutation (SHM) events across cells, and outputs an annotated phylogenetic tree showing how B cell clones evolved and expanded.
B cell lineage tracing is a tool for understanding how the immune system builds and refines antibody responses. By capturing the order and location of these mutations, we can get insights into the selection shaping BCR evolution and antibody creation.


## Background: ##
* During an immune response, B cells in our body undergo clonal expansion and affinity maturation in order to produce a large pool of high‑affinity, antigen‑specific antibodies that can neutralize the pathogens more effectively.
As part of this process, the cells accumulate somatic hypermutations (SHM) in their BCR sequences, specifically in the V(D)J region, to improve the antigen.

* B cell lineage trees map the evolutionary "family history" of B cell clones as they mutate and divide during an immune response - all originated from one B cell.
  
* The trees are essential for couple of reasons:

  -**Tracking Affinity Maturation:** Trees illustrate how B cells repeatedly mutate and select for stronger, more precise antibody binding against specific pathogens.

  -**Isotype Switching:** By tracking the genetic changes, the trees can reveal class-switch recombination (change in the antiboy type), showing how B cells change their functional properties over time.

  -**Therapeutic Antibody Discovery:** We can use the trees to trace mutated B cell sequences backwards, identifying the most potent, broadly neutralizing antibodies to isolate for therapeutic use.

  -**Disease & Vaccine Research:** Analyzing tree shapes reveals whether an immune response is generating new, evolving antibodies or just re-stimulating older, less effective memory cells.


## The tool: ##
### How is it working:
I developed this Python script to construct BCR clonal lineage trees based on the fundamental principle of clonal evolution. The tool maps the evolutionary "family history" of B cell clones by rooting each tree in its unmutated common ancestor (germline)—representing the cell’s original state before it encountered an antigen. As these B cells divide during an immune response, they accumulate unique genetic changes through the SHM process. By analyzing the hierarchy of shared mutations, the script reconstructs the branching order of the cells, where the length of each branch reflects the number of genetic changes acquired over time. To ensure high-fidelity results, the code incorporates critical biological constraints such as isotype switching (the irreversible transition between antibody types) and tissue connectivity, allowing for a professional visualization of how immune lineages mature and spread throughout different organs.

### Input:
The tool will take an **input of tabular data in .xlsx**  with one row per cell.
**Column order does not matter**. The tool matches columns by name (case-insensitive) and also recognises common aliases automatically.
#### *The columns:*
#### Required columns

| Column | Description |
|---|---|
| `clone_id` | Clone group identifier |
| `VDJ_sequence_H` | Heavy chain VDJ nucleotide sequence |

#### Strongly recommended

| Column | Description | What happens if absent |
|---|---|---|
| `germline_alignment_d_mask` | IMGT-masked germline reference sequence | A per-clone consensus is estimated from the observed sequences (less accurate) |

#### Optional — tree still builds; missing fields show as `Unknown`

| Column | Description |
|---|---|
| `cell_id` | Unique cell barcode (row index used as fallback) |
| `VDJ_sequence_L` | Light chain VDJ nucleotide sequence — when present, H+L are concatenated for higher phylogenetic accuracy |
| `c_call` | Isotype / constant-region call (e.g. `IGHM`, `IGHG1`) — used for node colours and CSR refinement |
| `sample_id` | Sample / tissue of origin — used for node colours in heavy-only datasets |
| `cluster_annotated` | Cell-type or GC-subset annotation — used for node shapes |
| `clone_count` | Number of cells in the clone |
| `VDJ_aa_sequence_H` | Heavy chain amino acid sequence |
| `VDJ_aa_sequence_L` | Light chain amino acid sequence |
| `v_call_h` / `d_call_h` / `j_call_h` | Heavy chain V/D/J gene calls |
| `v_call_l` / `j_call_l` | Light chain V/J gene calls |
| `mu_count_h` / `mu_count_l` | Somatic hypermutation counts |

* I decided to start with processed files rather than raw data because raw scBCR‑seq output typically comes in complex bioinformatics formats, and handling it can easily become a full project on its own. For this project, I prefered to focus on the downstream biology (lineage relationships, somatic hypermutation patterns, and class-switch recombination), rather than re-implementing well known upstream steps.
  
* The tool tries a list of recognised aliases for each role, so many non-standard column names are detected automatically.

### Output:
The outputs are written to a single folder (default: `bcr_lineage_output/`).

| File | Format | Description |
|---|---|---|
| `tree_<clone_id>.png` | PNG | Phylogenetic tree for each clone |
| `mutation_table.xlsx` | Excel | The mutation events, isotypes, cell type and more across all processed clones |

Example of the tree structure:

<img width="425" height="650" alt="image" src="https://github.com/user-attachments/assets/752dfbb7-862a-4431-87c5-cc24752ede08" />

### *The Tree:*

### *Mutation table:*

| Column | Description |
|---|---|
| `clone_id` | Clone this edge belongs to |
| `node` | Original cell ID or internal node name |
| `seq_label` | Short label as shown on the tree (`seq1`, `anc2`, …) |
| `parent` | Parent node name |
| `is_observed` | `True` for real cells, `False` for inferred ancestral nodes |
| `isotype` | Isotype of the child node |
| `sample_id` | Tissue of origin if the information is available |
| `cluster_annotated` | Cell-type annotation |
| `branch_length` | Hamming distance to parent |
| `number_nucleotides_changes` | Number of nucleotide changes on this edge |
| `nucleotides_changes` | Changes list e.g. `A42T;C81G` |
| `number_amino_acid_changes` | Number of amino acid changes |
| `amino_acid_changes` | Changes list e.g. `K14R` |







## The technicalities: ##
As for now, the requirements.txt will include:
* pandas
* openpyxl
* biopython
* ete3
* matplotlib
* scipy
* numpy
* tkinter (for GUI )

## Installation and running the project

**Requirements:** Python 3.9+

### 1. Install dependencies

pip install -r requirements.txt

### 2. Run the project

**Option A — Graphical Interface (recommended):**

python gui.py

This opens a window where you can upload your `.xlsx` / `.csv` / `.tsv` file,
select the clone to analyze, and view or export the results.

**Option B — Command Line:**

python run_lineage.py - input your_file.xlsx

### 3. Run the tests

pytest tests/


### Notes:
This project is part of the Python Programming Course at the Weizmann Institute of Science.

You can view the course main repository here: https://github.com/Code-Maven/wis-python-course-2026-03

### Sources I used to create the tree script:
* https://mixcr.com/mixcr/guides/b-cell-lineages-webinar/
* https://www.antibodysociety.org/wordpress/wp-content/uploads/2021/11/Immcantation-webinar-slides.pdf
* Weber, Leah L., et al. "Isotype-aware inference of B cell clonal lineage trees from single-cell sequencing data." Cell Genomics 4.9 (2024).



















---

## Algorithm overview

1. **Format detection** — if `VDJ_sequence_L` is present and non-empty the dataset is treated as *paired* (H+L concatenated); otherwise *heavy-only*.
2. **Germline (UCA) construction** — extracted from `germline_alignment_d_mask` if available; otherwise a per-clone positional consensus is computed from the observed sequences.
3. **Distance matrix** — pairwise alignment-based Hamming distance; handles unequal sequence lengths via global pairwise alignment.
4. **Neighbor-Joining tree** — built with BioPython's `DistanceTreeConstructor.nj()`; the germline is injected as the outgroup and becomes the tree root.
5. **Polytomy collapsing** — branches shorter than the collapse threshold (default `1e-6`) are merged into polytomies, reflecting the biological reality that short BCR sequences often cannot support full bifurcating resolution.
6. **Isotype-aware NNI refinement** *(optional, default on)* — a local nearest-neighbour interchange search penalises parent→child edges that violate the irreversibility of class-switch recombination (e.g. IgG → IgM is disallowed).
7. **Fitch parsimony reconstruction** — ancestral sequences are inferred at every internal node so the mutation table covers all edges, not just terminal branches.
8. **Visualisation and export** — tree PNG and mutation table written to the output directory.

---

## Project structure

```
Final Project new version/
│
├── run_gui.py                      ← launch the GUI (recommended)
├── run_cli.py                      ← command-line entry point
├── requirements.txt
│
├── bcr_lineage_tracer/             ← package (do not run files inside directly)
│   ├── constants.py                   Isotype CSR order, gap chars, marker cycle
│   ├── loader.py                      BCRTreeLoader — reads xlsx, detects format,
│   │                                  resolves columns, builds germlines
│   ├── tracer.py                      LineageTracer — NJ tree, Fitch parsimony,
│   │                                  polytomy collapse, isotype NNI
│   ├── visualization.py               plot_tree() — matplotlib cladogram
│   ├── pipeline.py                    run() — orchestrates all steps, writes outputs
│   ├── gui.py                         Tkinter graphical interface
│   ├── __main__.py                    CLI argument parser
│   └── __init__.py                    Public API re-exports
│
└── tests/
    ├── conftest.py                    Shared fixtures (no patient data required)
    ├── test_loader.py                 32 tests — column detection, format detection,
    │                                  germline construction, IMGT ungapping,
    │                                  error and warning handling
    ├── test_tracer.py                 24 tests — tree structure, polytomies,
    │                                  isotype NNI, Fitch reconstruction,
    │                                  mutation table, distances
    ├── test_pipeline.py               12 tests — end-to-end integration
    └── test_visualization.py          10 tests — matplotlib output, legend, PNG
```

---

## Installation

**Python 3.9 or higher required.**

### 1. Install dependencies

```bash
pip install -r requirements.txt
```

> **tkinter** (required for the GUI) ships with the official Python installer on Windows and macOS.
> On Linux you may need to install it separately:
> ```bash
> # Ubuntu / Debian
> sudo apt install python3-tk
>
> # Fedora
> sudo dnf install python3-tkinter
> ```

### 2. Run the project

**Option A — Graphical Interface (recommended)**

```bash
python run_gui.py
```

This opens the GUI where you can:
- Browse for your `.xlsx` file
- View how many clones were detected and how many are eligible (≥ 2 cells)
- Select a specific clone from a searchable scrollable list, or process all
- Adjust advanced options (see below)
- Watch the live log and open the output folder when done

**Option B — Command Line**

```bash
python run_cli.py --input your_file.xlsx
```

Common options:

| Flag | Default | Description |
|---|---|---|
| `--input` / `-i` | *(required)* | Path to `.xlsx` input file |
| `--output-dir` / `-o` | `bcr_lineage_output` | Folder for all output files |
| `--clone-id` | *(all clones)* | Process only one specific clone |
| `--max-clones N` | *(all)* | Process at most N clones — useful for quick tests |
| `--collapse-threshold` | `1e-6` | Branch-length cutoff for polytomy collapsing |
| `--no-isotype-refine` | *(refine on)* | Disable isotype-aware NNI |
| `--gui` | — | Launch the GUI instead |

### 3. Run the tests

```bash
pip install pytest
pytest tests/
```

All 78 tests run from synthetic in-memory data — no patient files needed.

---

## Advanced options explained

### Max clones *(default: all)*
Limits how many clonal families are processed in one run. Leave blank to process every clone. Useful for quick testing: set to `5` to verify the pipeline works on your file before committing to a full run.

### Collapse threshold *(default: `1e-6`)*
BCR sequences are short (~300–700 bp) and differ from each other by only a few mutations, so the tree algorithm often cannot distinguish whether two cells diverged before or after a mutation. Branches shorter than this threshold are collapsed into polytomies — which is biologically more honest than forcing an arbitrary ordering. The default removes only floating-point noise. Increase (e.g. `0.001`) for a simpler, more collapsed tree; decrease to `0` to keep every branch.

### Isotype-aware NNI refinement *(default: on)*
After building the initial tree from sequence distances, this step checks for parent→child edges that imply a biologically impossible isotype transition. Class-switch recombination (CSR) is irreversible — a cell that has switched to IgG cannot revert to IgM. The refinement performs local swaps to fix such violations, using the sequence signal as a tiebreaker so the topology changes as little as possible. Turn off if your file has no `c_call` column, or if you want the unmodified sequence-distance tree for comparison.

---

## Requirements

```
biopython>=1.80
matplotlib>=3.6
numpy>=1.23
openpyxl>=3.0
pandas>=1.5
```

---

## Biological background

This tool focuses on a stage of B cell biology that is **less commonly explored** in existing pipelines: the **within-clone evolutionary dynamics** after antigen encounter.

During a germinal centre (GC) reaction, B cells undergo rapid somatic hypermutation (SHM) and selection. Cells with beneficial mutations in their BCR expand; others are eliminated. Simultaneously, some cells undergo class-switch recombination (CSR), changing from IgM to IgG, IgA, or IgE to alter effector function.

By reconstructing the lineage tree of a clone — the branching order of cells and their accumulated mutations — this tool makes it possible to:

- Trace the **trajectory of affinity maturation** from the naïve germline to mature antibody-secreting plasma cells
- Identify **convergent evolution** (independent mutations at the same position in separate branches)
- Map **tissue migration** between organs (when `sample_id` encodes anatomical sites)
- Pinpoint the **branching point at which CSR occurred** in the lineage

The approach is grounded in methods from the Immcantation framework (Dowser / IgPhyML) and draws on concepts from TRIBAL (Tree Inference of B cell Clonal Lineages) for isotype-aware refinement.

*‏




