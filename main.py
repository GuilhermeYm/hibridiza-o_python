import random

import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem

moleculas = {
    "etano": "CC",
    "eteno": "C=C",
    "etino": "C#C",
    "propanona": "CC(=O)C",
    "propan-2-ol": "CC(O)C",
    "ácido etanoico": "CC(=O)O",
    "benzeno": "c1ccccc1",
    "propanal": "CCC=O",
    "acetonitrila": "CC#N",
}


def gerar_molecula_3d(smiles, nome_do_arquivo):
    print(nome_do_arquivo)
    if not smiles:
        raise ValueError("A string SMILES não pode ser vazia.")

    mol = Chem.MolFromSmiles(smiles)

    if not mol:
        raise ValueError("O valor inserido não é um SMILES válido.")

    mol = Chem.AddHs(mol)

    resultado_3d = AllChem.EmbedMolecule(mol)

    if resultado_3d != 0:
        raise ValueError("Não foi possível gerar a conformação 3D da molécula.")

    mol_block = Chem.MolToMolBlock(mol)

    view = py3Dmol.view(width=600, height=400)
    view.addModel(mol_block, "mol")
    view.setStyle({"stick": {}})
    view.zoomTo()

    html = view._make_html()

    nome_arquivo = nome_do_arquivo.removesuffix(".mol") + ".html"
    with open(nome_arquivo, "w", encoding="utf-8") as file:
        file.write(html)

    print("Arquivo criado com sucesso.")

    question = (
        input(
            "Você quer responder a um quiz sobre a hibridização dos átomos dessa molécula? (s/n) > "
        )
        .strip()
        .lower()
    )

    if question == "s":
        quiz_hibrização(smiles)


def quiz_hibrização(smiles):
    mol = Chem.MolFromSmiles(smiles)
    for atom in mol.GetAtoms():
        print(str(atom.GetHybridization()))
        if atom.GetSymbol() == "C":
            index_atom = atom.GetIdx()
            hybridization_atom = str(atom.GetHybridization()).strip().lower()

            resposta = (
                input(f"Qual é a hibridização desso carbono {index_atom}?")
                .strip()
                .lower()
            )
            if resposta == hybridization_atom:
                print("Parabéns!")
                continue
            else:
                print(f"Errado era {hybridization_atom}")


value = random.choice(list(moleculas.values()))
name = filter(lambda item: item[1] == value, moleculas.items())

gerar_molecula_3d(value, f"{list(name)[0][0]}.mol")
