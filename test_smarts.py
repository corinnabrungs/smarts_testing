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
        Chem.GetSSSR(self.smarts_mol)


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
                negative_smiles=["CC(=O)O", "CCC(=O)C"],
            )
        )
        test_structure_smarts(
            SmartsTest(
                description="aldehyde2",
                smarts="O=[CH1][CH2]C",
                positive_smiles=["CCC=O", "CCCCCC=O"],
                negative_smiles=["CC(=O)O", "CCC(=O)C"],
            )
        )
        test_structure_smarts(
            SmartsTest(
                description="terminal_alkene",
                smarts="[$([CX3H2]=[C]([#6,#1])[#6])]",
                positive_smiles=["CC(=C)C", "C(=C)C"],
                negative_smiles=["C=C", "CC=CC"],
            )
        )
        test_structure_smarts(
            SmartsTest(
                description="beta_lactone",
                smarts="O=C1OCC1",
                positive_smiles=["C1COC1=O", "CCCC1COC1=O"],
                negative_smiles=["CC1CCC(=O)O1"],
            )
        )
        test_structure_smarts(
            SmartsTest(
                description="enon",
                smarts="CC(=O)C([C,#1,O])=C([C,#1,O])[C,#1,O]",
                positive_smiles=["CC(=O)C=C", "CC(=O)C(O)=CC"],
                negative_smiles=["CCC=O", "C=CCC=O"],
            )
        )
        test_structure_smarts(
            SmartsTest(
                description="enacid",
                smarts="[C,#1]OC(=O)C([C,#1,O])=C([C,#1,O])[C,#1,O]",
                positive_smiles=["COC(=O)C=C(O)C"],
                negative_smiles=["CCC=O"],
            )
        )
        test_structure_smarts(
            SmartsTest(
                description="epoxyketone",
                smarts="O=C([#6])C1C([#6])O1",
                positive_smiles=["O=C(C)C1C(C)O1", "O=C(c1ccccc1)C1C(C)O1"],
                negative_smiles=["CC=O"],
            )
        )
        test_structure_smarts(
            SmartsTest(
                description="anthraquinone",
                smarts="c12ccccc1C(=O)c3ccccc3C2=O",
                positive_smiles=[
                    "c12cc(C)c(O)c(C)c1C(=O)c3c(C)cccc3C2=O",
                    "COC1=CC=CC2=C1C(=O)C3=C(C2=O)C=CC=C3O",
                    "CC1=CC2=C(C(=C1)O)C(=O)C3=C(C2=O)C=C(C=C3O)O",
                ],
                negative_smiles=["CC=O"],
            )
        )
        test_structure_smarts(
            SmartsTest(
                description="naphtochinon",
                smarts="O=C-1-C([C,#1,O])=C([C,#1,O])-C(=[O])-c:2:c:c:c:c:c-1:2",
                positive_smiles=[
                    "O=C-1-C(O)=C(C)-C(=O)-c:2:c:c:c:c:c-1:2",
                    "C1=CC=C2C(=C1)C(=O)C=C(C2=O)C",
                ],
                negative_smiles=["CC=O"],
            )
        )
        test_structure_smarts(
            SmartsTest(
                description="prim_sec_amine",
                smarts="[C]N([#1])[C,#1]",
                positive_smiles=[
                    "CN",
                    "CNC",
                ],
                negative_smiles=["CC=O"],
            )
        )
        test_structure_smarts(
            SmartsTest(
                description="hydroxy-amide",
                smarts="CC(=O)N([OH])C",
                positive_smiles=[
                    "CC(=O)N(O)C",
                ],
                negative_smiles=["CC=O"],
            )
        )
        test_structure_smarts(
            SmartsTest(
                description="n_hydroxy_pyrrol",
                smarts="[OH]n1cccc1",
                positive_smiles=[
                    "On1cccc1",
                ],
                negative_smiles=["CC=O"],
            )
        )
        print("All smarts tested successfully")
