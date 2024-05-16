### Testing of SMARTS substructures by positive and negative smiles list

Adding new tests is simple. Open test_smarts.py and add a new SMARTS to test like below. Then run unit

```python
test_structure_smarts(
    SmartsTest(
        description="terminal_alkene",
        smarts="[$([CX3H2]=[C]([#6,#1])[#6])]",
        # List of structures that contain the SMARTS
        positive_smiles=["CC(=C)C", "C(=C)C"],
        # List of structures that should NOT match the SMARTS
        negative_smiles=["C=C", "CC=CC"],
    )
)
```

### Run tests

```bash
python -m unittest
```


### Example error

The test will print that all SMARTS were tested successfully if all tests are correct. The first failing Test will exit the program with a message like this, showing the failing SMILES and SMARTS pair and expected result: match/no match:

```
Traceback (most recent call last):
  File "C:\git\smarts_testing\test_smarts.py", line 62, in test_smarts
    test_structure_smarts(
  File "C:\git\smarts_testing\test_smarts.py", line 29, in test_structure_smarts
    raise ValueError(
ValueError: terminal_alkene, had match for negative smiles test: CC=C and smarts [$([CX3H2]=[C]([#6,#1])[#6])]

```
