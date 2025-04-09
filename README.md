# Energy Analysis: Spike RBD–ACE2 Protein-Protein Interface

## 📌 Description

This project explores the interface of the Spike RBD–ACE2 protein complex. Our goal is to understand how individual residues contribute to binding energy—one of the main factors explaining the interaction stability within this complex.

---

## 🚀 Installation and Setup

To replicate our pipeline and obtain the same results, follow the steps outlined below:

### Preparation: Cleaning the Structure

1. **Download** the PDB structure for `6m0j`.
2. **Install dependencies**:
   ```bash
   pip install biobb_structure_checking biopython
   ```
3. **Run** the cleaning script:
   ```bash
   python preparation.py
   ```
   - Finds chains
   - Removes heteroatoms
   - Performs a quality check

---

## 🔬 Step-by-Step Pipeline

### Step 1: Finding Interface Residues

- Open `6m0j_fixed` in **PyMOL**.
- Identify the largest distance between potentially interacting residues across chains.
- Add 2 Å to the measured distance to include adjacent residues.
- Run:
  ```bash
  python step1.py
  ```
  - Outputs a list of interface residues in both chains.
- Reference video: [PyMOL interface residue tutorial](https://www.youtube.com/watch?v=hcnnKrlqa9M)

---

### Step 2: Finding Interaction Energies

- **Download and install**:
  - `vdwprm`, `aaLib.lib` (from GitHub/Moodle)
  - `NACCESS`
  - Install NACCESS:
    ```bash
    csh install.scr
    ```
- Run:
  ```bash
  python step2.py
  ```
  - Calculates:
    - Total complex energy
    - Interface energy
  - Compares energies and extracts conclusions

---

### Step 3: Alanine Scanning Experiment

- Run:
  ```bash
  python step3.py
  ```
  - Mutates interface residues to alanine (ALA)
  - Calculates energy changes for each mutation
  - Outputs a table of energy differences
  - Generates a barplot for visualization

---

### 🖼️ Step 4: Generating PyMOL Images

- Open `6m0j_fixed` in **PyMOL**
- Follow the provided pipeline (`pymol_pipeline.txt`) to:
  - Highlight interface
  - Capture images of mutations and bond changes
- Tutorial video: [PyMOL visualization guide](https://www.youtube.com/watch?v=wGRMGEnHPdg)

---

## 👥 Contributors
 
- María López  [MLM19](https://github.com/MLM19)
- Stephanie Padilla [MstephPad](https://github.com/MstephPad)
