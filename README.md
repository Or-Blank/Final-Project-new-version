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

#### *The Tree:*
- Nodes are labelled **seq1, seq2, seq3 …** (observed cells, top-to-bottom) and **anc1, anc2 …** (inferred ancestral nodes).
- Node **colour** encodes isotype (`c_call`) for paired H+L datasets, or tissue of origin (`sample_id`) for heavy-only datasets.
- Node **shape** encodes cell-type annotation (`cluster_annotated`).
- The **germline root** is always drawn as a large black square labelled "Germline".
- The x-axis shows **cumulative mutation distance from the germline**.

#### *The Mutation table:*

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

* Note: There is a recurring pattern in both the tree and the table where several cell-type annotations appear to represent essentially the same population (e.g., 2‑B cell, 3‑B cell). These are all GC B cells, but they differ slightly in characteristics that are not reflected in the tree or the table and become apparent only in other experiments.


## The technicalities: ##
The requirements.txt include:
* biopython
* matplotlib
* numpy
* openpyxl
* pandas
* pytest
* tkinter (for GUI )

## Installation and running the project

**Requirements:** Python 3.9+

### 1. Install dependencies

```bash
pip install -r requirements.txt
```

### 2. Run the project
```bash
python run_gui.py
```
This opens a window where you can upload your `.xlsx` file,
select the clone to analyze, and view the results.

### 3. Run the tests (optional)
```bash
pytest tests/
```


### Notes:
This project is part of the Python Programming Course at the Weizmann Institute of Science.

You can view the course main repository here: https://github.com/Code-Maven/wis-python-course-2026-03


### Sources I used to create the tree script:
* https://mixcr.com/mixcr/guides/b-cell-lineages-webinar/
* https://www.antibodysociety.org/wordpress/wp-content/uploads/2021/11/Immcantation-webinar-slides.pdf
* Weber, Leah L., et al. "Isotype-aware inference of B cell clonal lineage trees from single-cell sequencing data." Cell Genomics 4.9 (2024).
