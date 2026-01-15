# Contributing to MolecularDFT

First off, thank you for considering contributing to MolecularDFT! 🎉

This document provides guidelines and steps for contributing. Following these guidelines helps communicate that you respect the time of the developers managing this project.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [How Can I Contribute?](#how-can-i-contribute)
- [Development Setup](#development-setup)
- [Style Guidelines](#style-guidelines)
- [Pull Request Process](#pull-request-process)

---

## Code of Conduct

This project and everyone participating in it is governed by our commitment to providing a welcoming and inclusive environment. By participating, you are expected to uphold this standard. Please be respectful and constructive in all interactions.

---

## Getting Started

### Issues

- **Bug Reports**: Use the bug report template to describe what went wrong
- **Feature Requests**: Use the feature request template to propose new functionality
- **Questions**: Open a discussion for general questions

Before creating an issue, please:
1. Search existing issues to avoid duplicates
2. Use a clear, descriptive title
3. Provide as much relevant information as possible

---

## How Can I Contribute?

### 🐛 Reporting Bugs

When filing a bug report, please include:

```markdown
**Environment:**
- OS: [e.g., Ubuntu 22.04]
- Python version: [e.g., 3.10.6]
- Package version: [e.g., 1.0.0]

**Description:**
A clear description of the bug.

**Steps to Reproduce:**
1. Step one
2. Step two
3. ...

**Expected Behavior:**
What you expected to happen.

**Actual Behavior:**
What actually happened.

**Screenshots/Logs:**
If applicable, add screenshots or error logs.
```

### 💡 Suggesting Features

Feature requests are welcome! Please include:

- **Use case**: Why is this feature needed?
- **Proposed solution**: How should it work?
- **Alternatives considered**: Other approaches you've thought about
- **Additional context**: Any other relevant information

### 🔧 Code Contributions

We welcome contributions in:

- **Bug fixes**
- **New features**
- **Documentation improvements**
- **Test coverage**
- **Performance optimizations**

---

## Development Setup

### 1. Fork and Clone

```bash
# Fork via GitHub UI, then:
git clone https://github.com/YOUR_USERNAME/MolecularDFT.git
cd MolecularDFT
```

### 2. Create Development Environment

```bash
# Using conda (recommended)
conda create -n moleculardft-dev python=3.10
conda activate moleculardft-dev

# Install dependencies
conda install -c conda-forge rdkit pyscf numpy scipy pandas matplotlib

# Install in development mode
pip install -e ".[dev]"
```

### 3. Create a Branch

```bash
# For features
git checkout -b feature/your-feature-name

# For bug fixes
git checkout -b fix/issue-description

# For documentation
git checkout -b docs/what-you-changed
```

### 4. Make Your Changes

- Write clean, documented code
- Add tests for new functionality
- Update documentation as needed

### 5. Run Tests

```bash
# Run all tests
pytest tests/

# Run with coverage
pytest --cov=src tests/

# Run specific test file
pytest tests/test_molecule.py
```

---

## Style Guidelines

### Python Code Style

We follow [PEP 8](https://peps.python.org/pep-0008/) with these specifics:

```python
# Good: Descriptive names
def calculate_molecular_descriptors(molecule: DrugMolecule) -> Dict[str, float]:
    """
    Calculate comprehensive molecular descriptors.
    
    Parameters
    ----------
    molecule : DrugMolecule
        The molecule to analyze.
        
    Returns
    -------
    Dict[str, float]
        Dictionary of descriptor names to values.
        
    Examples
    --------
    >>> mol = DrugMolecule("aspirin", smiles="CC(=O)OC1=CC=CC=C1C(=O)O")
    >>> descriptors = calculate_molecular_descriptors(mol)
    >>> print(descriptors['MW'])
    180.16
    """
    pass
```

### Commit Messages

Use conventional commits:

```
type(scope): subject

body (optional)

footer (optional)
```

Types:
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation
- `style`: Formatting (no code change)
- `refactor`: Code restructuring
- `test`: Adding tests
- `chore`: Maintenance

Examples:
```
feat(dft): add D4 dispersion correction support

fix(admet): correct Lipinski HBD calculation

docs(readme): update installation instructions
```

### Documentation

- Use NumPy-style docstrings
- Include type hints
- Add examples where helpful

---

## Pull Request Process

### 1. Before Submitting

- [ ] Code follows style guidelines
- [ ] Tests pass locally
- [ ] Documentation is updated
- [ ] Commit messages follow conventions
- [ ] Branch is up to date with `main`

### 2. PR Description Template

```markdown
## Description
Brief description of changes.

## Type of Change
- [ ] Bug fix (non-breaking change fixing an issue)
- [ ] New feature (non-breaking change adding functionality)
- [ ] Breaking change (fix or feature causing existing functionality to change)
- [ ] Documentation update

## Testing
Describe tests performed.

## Checklist
- [ ] My code follows the project style guidelines
- [ ] I have performed a self-review
- [ ] I have commented my code where necessary
- [ ] I have updated the documentation
- [ ] My changes generate no new warnings
- [ ] I have added tests that prove my fix/feature works
- [ ] New and existing tests pass locally
```

### 3. Review Process

1. Submit PR against `main` branch
2. Automated checks will run
3. Maintainer will review within 1-2 weeks
4. Address any feedback
5. Once approved, maintainer will merge

---

## Recognition

Contributors will be recognized in:
- The project README
- Release notes
- CONTRIBUTORS.md file

---

## Questions?

Feel free to open a discussion or reach out to the maintainers.

Thank you for contributing! 🚀
