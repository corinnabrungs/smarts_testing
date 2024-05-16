# Testing of SMARTS substructures by positive and negative smiles list

Adding new tests is simple. Open test_smarts.py and add a new SMARTS to test like below:

```python
test_structure_smarts(
    SmartsTest(
        description="terminal_alkene",
        smarts="[$([CX3H2]=[C]([#6,#1])[#6])]",
        positive_smiles=["CC(=C)C", "C(=C)C"],  # List of structures that contain the SMARTS
        negative_smiles=["C=C", "CC=CC"],                # List of structures that should not match the SMARTS
    )
)
```
