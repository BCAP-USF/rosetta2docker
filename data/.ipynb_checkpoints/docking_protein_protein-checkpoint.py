from pyrosetta import *
from pyrosetta.rosetta.protocols.docking import *
from pyrosetta.rosetta.core.scoring import *

# 1. Inicializar o PyRosetta
init()

# 2. Carregar e combinar os dois PDBs individuais
pose_receptor = pose_from_file("RAR.pdb")
pose_ligante = pose_from_file("Curcumina.sdf")

pose_inicial = Pose()
pose_inicial.assign(pose_receptor)
pose_inicial.append_pose_to_pose(pose_ligante, jump_new_ider=True)

# 3. Configurar a FoldTree para o movimento (Cadeia A fixa, Cadeia B móvel)
setup_foldtree(pose_inicial, "A_B", Vector1())

# 4. Definir as funções de pontuação (ScoreFunctions)
scorefxn_high = create_score_function("ref2015")
scorefxn_low = create_score_function("interchain_cen")

# 5. Instanciar os motores de movimentação (Movers)
to_centroid = SwitchResidueTypeSetMover("centroid")
to_fullatom = SwitchResidueTypeSetMover("fa_standard")
dock_lowres = DockingLowRes(scorefxn_low, 1)   # 1 representa o índice do jump
dock_highres = DockingHighRes(scorefxn_high, 1)

# 6. Configurar o Distribuidor de Trabalhos em Lote (PyJobDistributor)
# Parâmetros: Nome do prefixo dos arquivos, total de decoys, nome do arquivo de scores
nome_prefixo = "docking_decoy"
total_decoys = 100
score_file = "docking_scores.fac"

jd = PyJobDistributor(nome_prefixo, total_decoys, score_file)

print(f"Iniciando a produção em lote de {total_decoys} decoys...")

# 7. Loop de execução gerenciado pelo Job Distributor
while not jd.job_done():
    # Cria uma cópia limpa da estrutura inicial para cada nova tentativa
    test_pose = pose_inicial.clone()
    
    # Fase de Baixa Resolução (Centroid)
    to_centroid.apply(test_pose)
    dock_lowres.apply(test_pose)
    
    # Retorno para Alta Resolução (Full-Atom) e Refinamento
    to_fullatom.apply(test_pose)
    dock_highres.apply(test_pose)
    
    # O Job Distributor avalia o score final e exporta o PDB automaticamente
    jd.output_decoy(test_pose)

print("Processamento em lote concluído com sucesso!")
