import pyrosetta
pyrosetta.init()

from pyrosetta.rosetta.core.chemical import ResidueType, sdf, write_topology_file


def gerar_params(sdf_path, nome="LG1", params_path="LIG.params"):
    """Gera o .params a partir de um .sdf usando PyRosetta."""
    for mutable_type in sdf.convert_to_ResidueTypes(sdf_path):
        mutable_type.name(nome)
        # O PyRosetta deriva name3 do título da molécula ("969516" -> "969"),
        # e é por ele que write_topology_file monta o IO_STRING. Sem isto o
        # .params declara um resíduo que não bate com o PDB e o docking falha
        # com "Unrecognized residue".
        mutable_type.name3(nome)
        mutable_type.name1(nome[-1])
        residue_type = ResidueType.make(mutable_type)
        write_topology_file(residue_type, params_path)
        print(f"Arquivo {params_path} gerado com sucesso!")
        return params_path
    raise RuntimeError(f"Nenhuma molécula válida encontrada em {sdf_path}")


def gerar_pdb(sdf_path, nome="LG1", pdb_path="LIG_0001.pdb", params_path=None):
    """Gera o PDB do ligante a partir do bloco de átomos do .sdf.

    O write_topology_file do PyRosetta só emite topologia, então o PDB de
    coordenadas é montado lendo o molfile diretamente. Os nomes de átomo
    vêm do .params para garantir que PDB e topologia usem os mesmos.
    """
    linhas = open(sdf_path).read().splitlines()
    n_atoms = int(linhas[3][:3])
    bloco = linhas[4:4 + n_atoms]

    if params_path is None:
        raise RuntimeError("Informe params_path: os nomes de átomo vêm dele.")
    nomes = [l.split()[1] for l in open(params_path) if l.startswith("ATOM ")]
    if len(nomes) != n_atoms:
        raise RuntimeError(
            f"{n_atoms} átomos no SDF mas {len(nomes)} no .params — "
            f"PDB e topologia ficariam inconsistentes.")

    with open(pdb_path, "w") as f:
        for i, (linha, nome_atom) in enumerate(zip(bloco, nomes), start=1):
            x, y, z, elem = linha.split()[:4]
            f.write(
                f"HETATM{i:5d} {nome_atom:<4s} {nome:<3s} X{i:4d}    "
                f"{float(x):8.3f}{float(y):8.3f}{float(z):8.3f}"
                f"  1.00  0.00          {elem.upper():<2s}\n"
            )
        f.write("TER\nEND\n")
    print(f"Arquivo {pdb_path} gerado com sucesso!")
    return pdb_path


def gerar_sdf(pdb_path, sdf_origem, sdf_path):
    """Gera o .sdf do ligante a partir do PDB gerado por gerar_pdb.

    O PDB não guarda ordens de ligação (nem gerar_pdb escreve CONECT), então
    o bloco de bonds — e todos os demais campos do molfile — vêm do SDF de
    origem. Como gerar_pdb escreve os átomos na mesma ordem do SDF, basta
    substituir as coordenadas de cada linha do bloco de átomos; o resto é
    copiado sem alteração. Os elementos são conferidos átomo a átomo para
    falhar alto se as ordens divergirem.
    """
    linhas = open(sdf_origem).read().splitlines()
    n_atoms = int(linhas[3][:3])
    bloco = linhas[4:4 + n_atoms]
    if len(bloco) != n_atoms:
        raise RuntimeError(
            f"SDF de origem truncado: {n_atoms} átomos declarados mas "
            f"{len(bloco)} linhas no bloco de átomos.")

    atoms_pdb = [l for l in open(pdb_path) if l.startswith("HETATM")]
    if len(atoms_pdb) != n_atoms:
        raise RuntimeError(
            f"{len(atoms_pdb)} átomos no PDB mas {n_atoms} no SDF de origem — "
            f"as coordenadas não poderiam ser alinhadas pela ordem.")

    for i, (linha_pdb, linha_sdf) in enumerate(zip(atoms_pdb, bloco), start=1):
        # Colunas fixas dos dois formatos — nunca split() em linha de PDB.
        elem_pdb = linha_pdb[76:78].strip().upper()
        elem_sdf = linha_sdf.split()[3].upper()
        if elem_pdb != elem_sdf:
            raise RuntimeError(
                f"Átomo {i}: elemento {elem_pdb or '?'} no PDB mas {elem_sdf} "
                f"no SDF de origem — a ordem dos átomos divergiu.")

    with open(sdf_path, "w") as f:
        # título, programa, comentário e counts do SDF de origem
        f.writelines(f"{l}\n" for l in linhas[:4])
        for linha_pdb, linha_sdf in zip(atoms_pdb, bloco):
            x = float(linha_pdb[30:38])
            y = float(linha_pdb[38:46])
            z = float(linha_pdb[46:54])
            # Substitui só os 3 campos de coordenada (%10.4f); símbolo e
            # demais campos V2000 da linha original são preservados.
            f.write(f"{x:10.4f}{y:10.4f}{z:10.4f}{linha_sdf[30:]}\n")
        # bonds, M  END, propriedades e $$$$ do SDF de origem
        f.writelines(f"{l}\n" for l in linhas[4 + n_atoms:])
    print(f"Arquivo {sdf_path} gerado com sucesso!")
    return sdf_path


if __name__ == "__main__":
    import sys
    gerar_params(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "LG1")
