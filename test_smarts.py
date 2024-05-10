from unittest import TestCase
import rdkit.Chem as Chem

from dataclasses import dataclass, field


@dataclass
class SmartsTest:
    description: str
    smarts: str
    positive_smiles: list[str]
    negative_smiles: list[str]
    smarts_mol: Chem.rdchem.Mol = field(init=False)

    def __post_init__(self):
        self.smarts_mol = Chem.MolFromSmarts(self.smarts)


def test_structure_smarts(case: SmartsTest):
    for positive in case.positive_smiles:
        if count_substructure_matches(case.smarts_mol, positive) == 0:
            raise ValueError(
                f"{case.description}, had NO match for positive smiles test: {positive} and smarts {case.smarts}"
            )

    for negative in case.negative_smiles:
        if count_substructure_matches(case.smarts_mol, negative) != 0:
            raise ValueError(
                f"{case.description}, had match for negative smiles test: {negative} and smarts {case.smarts}"
            )


def count_substructure_matches(smarts_mol, smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
    except Exception as e:
        print(f"could not convert smiles to mol: {smiles}")
        raise e
    return len(mol.GetSubstructMatches(smarts_mol))


class Test(TestCase):
    def test_smarts(self):
        test_structure_smarts(
            SmartsTest(
                description="aldehyde",
                smarts="[CX3H1](=O)[CH2]C",
                positive_smiles=["CCC=O"],
                negative_smiles=["CC(=O)OH", "CC(=O)C"],
            )
        )
