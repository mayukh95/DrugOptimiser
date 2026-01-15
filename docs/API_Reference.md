# MolecularDFT API Reference

## Table of Contents

- [Core Module](#core-module)
- [Quantum Module](#quantum-module)
- [Interactive Module](#interactive-module)
- [Analysis Module](#analysis-module)
- [Visualization Module](#visualization-module)
- [Utils Module](#utils-module)

---

## Core Module

### `DrugMolecule`

The fundamental molecular data structure with RDKit integration.

```python
from src.core import DrugMolecule

mol = DrugMolecule(
    name="Aspirin",
    smiles="CC(=O)OC1=CC=CC=C1C(=O)O",
    cas="50-78-2"  # optional
)
```

#### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | `str` | Molecule identifier |
| `smiles` | `str` | SMILES string |
| `mol` | `rdkit.Chem.Mol` | RDKit molecule object |
| `elements` | `List[str]` | Atomic symbols |
| `initial_coords` | `np.ndarray` | 3D coordinates (N×3) |
| `properties` | `Dict[str, float]` | Calculated descriptors |

#### Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `calculate_descriptors()` | `Dict` | Calculate 50+ molecular descriptors |
| `to_xyz()` | `str` | Export as XYZ format |
| `to_sdf()` | `str` | Export as SDF format |

---

### `ADMETPredictor`

Comprehensive ADMET property prediction.

```python
from src.core import ADMETPredictor

predictor = ADMETPredictor()
results = predictor.predict_all(mol)
```

#### Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `check_lipinski(mol)` | `Dict` | Lipinski Rule of 5 |
| `check_veber(mol)` | `Dict` | Veber's rules |
| `check_ghose(mol)` | `Dict` | Ghose filter |
| `check_egan(mol)` | `Dict` | Egan filter |
| `check_muegge(mol)` | `Dict` | Muegge filter |
| `check_pains(mol)` | `Dict` | PAINS alerts |
| `check_brenk(mol)` | `Dict` | Brenk alerts |
| `calculate_sa_score(mol)` | `float` | Synthetic Accessibility |
| `predict_all(mol)` | `Dict` | All predictions |

---

## Quantum Module

### `EnhancedDFTCalculator`

PySCF-based DFT calculator with extensive options.

```python
from src.quantum import EnhancedDFTCalculator

calculator = EnhancedDFTCalculator(
    functional='B3LYP',
    basis='6-31G*',
    dispersion='d3bj',
    solvation='pcm',
    solvent='water',
    max_scf_cycles=100
)
```

#### Constructor Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `functional` | `str` | `'B3LYP'` | DFT functional |
| `basis` | `str` | `'6-31G*'` | Basis set |
| `dispersion` | `str` | `None` | D3, D3BJ, D4 |
| `solvation` | `str` | `None` | PCM, COSMO, ddCOSMO |
| `solvent` | `str` | `'water'` | Solvent name |
| `max_scf_cycles` | `int` | `100` | Max SCF iterations |
| `scf_damping` | `float` | `0.0` | SCF damping factor |
| `level_shift` | `float` | `0.0` | Level shift (Ha) |

#### Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `calculate_energy(coords, elements, charge, mult)` | `float` | Single-point energy |
| `calculate_gradient(coords, elements, charge, mult)` | `np.ndarray` | Energy gradient |
| `save_wavefunction(...)` | `bool` | Save NPZ + cube files |

---

### `ClassicalGeometryOptimizer`

Geometry optimization using scipy algorithms.

```python
from src.quantum import ClassicalGeometryOptimizer

optimizer = ClassicalGeometryOptimizer(
    dft_calculator=calculator,
    method='L-BFGS-B'
)

optimized_coords, history = optimizer.optimize(
    initial_coords=mol.initial_coords,
    elements=mol.elements,
    charge=0,
    multiplicity=1,
    max_iter=50,
    force_tol=0.001
)
```

#### Supported Methods
- `BFGS`, `L-BFGS-B`, `CG`, `Newton-CG`

---

## Interactive Module

### `DFTControlPanel`

Single-molecule DFT optimization interface.

```python
from src.interactive import DFTControlPanel

panel = DFTControlPanel()
panel.load_molecules(molecules)
display(panel.display())
```

### `BatchDFTControlPanel`

Parallel batch optimization interface.

```python
from src.interactive import BatchDFTControlPanel

batch_panel = BatchDFTControlPanel()
batch_panel.load_molecules(molecules)
display(batch_panel.display())
```

### `InteractivePropertyFilter`

ADMET property filtering with histograms.

```python
from src.interactive import InteractivePropertyFilter

filter_widget = InteractivePropertyFilter(molecules)
display(filter_widget.display())

# Get filtered molecules
filtered = filter_widget.get_filtered_molecules()
```

### `MoleculeComparator`

Side-by-side molecule comparison.

```python
from src.interactive import MoleculeComparator

comparator = MoleculeComparator(molecules)
display(comparator.display())
```

---

## Analysis Module

### `DFTDataCollector`

Aggregate results from wavefunction files.

```python
from src.analysis import DFTDataCollector

collector = DFTDataCollector(output_dir='./optimized_molecules')
df = collector.collect_data()
df.to_csv('results.csv', index=False)
```

#### Collected Properties

| Property | Description |
|----------|-------------|
| `Total_Energy_Ha` | Total electronic energy |
| `HOMO_eV` | HOMO orbital energy |
| `LUMO_eV` | LUMO orbital energy |
| `Gap_eV` | HOMO-LUMO gap |
| `Hardness_eV` | Chemical hardness |
| `Electrophilicity_Index` | Electrophilicity |
| `Max_Fukui_f_plus` | Max nucleophilic susceptibility |
| `Max_Spin_Density` | Radical character |
| `PSA_A2` | Polar surface area |
| `ESP_Min_au` | Minimum ESP |
| `ESP_Max_au` | Maximum ESP |

### `OrbitalVisualizer`

Molecular orbital visualization.

```python
from src.analysis import OrbitalVisualizer

viz = OrbitalVisualizer()
viz.load_wavefunction('path/to/wavefunction.npz')
fig = viz.plot_mo_diagram()
```

---

## Visualization Module

### `MolecularViewer3D`

3D molecular structure visualization.

```python
from src.visualization import MolecularViewer3D

viewer = MolecularViewer3D()
fig = viewer.plot_molecule(mol)
fig.show()
```

### `PropertyPlotter`

Property distribution plots.

```python
from src.visualization import PropertyPlotter

plotter = PropertyPlotter()
fig = plotter.plot_descriptor_distributions(molecules)
```

---

## Utils Module

### `export_molecules_to_sdf`

Export molecules to SDF format.

```python
from src.utils import export_molecules_to_sdf

export_molecules_to_sdf(molecules, 'output.sdf')
```

### `load_xyz_trajectory`

Load XYZ trajectory files.

```python
from src.utils import load_xyz_trajectory

frames = load_xyz_trajectory('trajectory.xyz')
```

---

## Configuration

### `Config`

Global configuration object.

```python
from src.core import Config

config = Config()
config.DEFAULT_FUNCTIONAL = 'PBE0'
config.DEFAULT_BASIS = 'def2-TZVP'
config.MAX_WORKERS = 4
```

---

## Error Handling

All modules raise descriptive exceptions:

```python
from src.core.exceptions import (
    MoleculeError,
    CalculationError,
    ConvergenceError
)

try:
    calculator.calculate_energy(coords, elements)
except ConvergenceError as e:
    print(f"SCF did not converge: {e}")
```
