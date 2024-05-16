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
    def test_extract_smarts_from_smiles(self):
        from rdkit.Chem import rdFMCS

        mols = [Chem.MolFromSmiles("CC(=C)C"), Chem.MolFromSmiles("CC(=C)")]
        smarts1 = rdFMCS.FindMCS(mols).smartsString
        smarts2 = rdFMCS.FindMCS(mols, ringMatchesRingOnly=True).smartsString
        print(f"{smarts1}+{smarts2}")

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
                negative_smiles=["C=C"],
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
                smarts="CC(C=C)=O",
                positive_smiles=["CC(=O)C=C"],
                negative_smiles=["CCC"],
            )
        )
        # still todo
        # test_structure_smarts(
        #     SmartsTest(
        #         description="enacid",
        #         smarts="O=C(O[C,H])C([C,H,O])=C([C,H,O])[C,H,O]",
        #         positive_smiles=[""],
        #         negative_smiles=[""],
        #     )
        # )
        # test_structure_smarts(
        #     SmartsTest(
        #         description="epoxyketone",
        #         smarts="O=C([#6])C1C([#6])O1",
        #         positive_smiles=[""],
        #         negative_smiles=[""],
        #     )
        # )
        # test_structure_smarts(
        #     SmartsTest(
        #         description="anthraquinone",
        #         smarts="c12c([*])c([*])c([*])c([*])c1C(=O)c3c([*])c([*])c([*])c([*])c3C2=O",
        #         smarts="c12ccccc1C(=O)c3ccccc3C2=O"
        #         positive_smiles=[""],
        #         negative_smiles=[""],
        #     )
        # )
        # test_structure_smarts(
        #     SmartsTest(
        #         description="naphtochinon",
        #         smarts="O=C-1-C([C,H,O])=C([C,H,O])-C(=[O])-c:2:c([*]):c([*]):c([*]):c([*]):c-1:2",
        #         positive_smiles=[""],
        #         negative_smiles=[""],
        #     )
        # )
