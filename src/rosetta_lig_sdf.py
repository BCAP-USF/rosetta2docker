#!/usr/bin/env python3
# =====================================================================
#  PIPELINE DE DOCKING PROTEÍNA-LIGANTE COM ROSETTA (RosettaLigand)
#  Plataforma-alvo: macOS / Apple Silicon (M1)
#
#  MODOS DE ANÁLISES
#  -----------------------------------------------
#
#   MODO "GLOBAL"  -> blind docking por MULTI-START. Detecta cavidades
#                     candidatas por heurística de enterramento, translada
#                     o ligante para cada centro e roda uma corrida
#                     INDEPENDENTE em cada um. Depois compara todas.
#                     É a forma canônica de busca cega no Rosetta.
#
#   MODO "LOCAL"   -> você informa as coordenadas do centro do sítio
#                     (CENTRO_SITIO). Uma única corrida, mais poses,
#                     refinamento fino. Use quando já souber o bolsão.
#
#  Em ambos os modos o ligante é EXPLICITAMENTE transladado para o ponto
#  de partida antes de montar o complexo.
# =====================================================================

import os
import re
import csv
import math
import random
import shutil
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

from molfile_to_params import gerar_params, gerar_pdb

# =====================================================================
# --------------------- 1. PAINEL DE CONTROLE -------------------------
# =====================================================================

# ---- Estratégia de busca --------------------------------------------
MODO_BUSCA = "GLOBAL"      # "GLOBAL" = blind docking multi-start
                           # "LOCAL"  = sítio conhecido (preencha CENTRO_SITIO)

# Coordenadas do centro do bolsão, em angstroms, no referencial do PDB
# do receptor. Só é lido quando MODO_BUSCA == "LOCAL".
# Como obter: centro de massa de um ligante cristalográfico do mesmo
# alvo, saída do fpocket/CASTp, ou centroide dos resíduos catalíticos.
CENTRO_SITIO = (0.0, 0.0, 0.0)

# ---- Orçamento de amostragem ----------------------------------------
# No modo GLOBAL o custo total é N_PONTOS_PARTIDA x POSES_POR_PONTO por
# par proteína-ligante. Estratégia de dois estágios recomendada:
#   triagem  -> N_PONTOS_PARTIDA=12, POSES_POR_PONTO=25  (300 poses)
#   refino   -> MODO_BUSCA="LOCAL" no melhor ponto, POSES_LOCAL=500
N_PONTOS_PARTIDA = 5       # quantas cavidades candidatas testar (modo GLOBAL)
POSES_POR_PONTO  = 5       # poses por ponto de partida (modo GLOBAL)
POSES_LOCAL      = 300     # poses na corrida única (modo LOCAL)

# ---- Recursos da máquina --------------------------------------------
# ATENÇÃO (MacBook Air M1): a máquina é fanless. Com 4 núcleos saturados
# por horas há throttling térmico severo e o tempo real fica muito acima
# da estimativa linear. Em corridas longas prefira 3.
MAX_NUCLEOS = 4

# ---- Preparação do receptor -----------------------------------------
# HETATM que devem ser MANTIDOS no receptor (íons metálicos, cofatores
# estruturais). Tudo que não estiver aqui é descartado, inclusive águas.
# Se o seu sítio tem Zn2+ catalítico e você esquecer de listá-lo aqui,
# o docking é fisicamente sem sentido.
HETATM_MANTER = {"ZN", "MG", "CA", "MN", "FE", "CU", "NA", "K"}

# Manter apenas o primeiro altloc de cada átomo (conformações
# alternativas do cristal geram átomos duplicados que o Rosetta trata mal)
FILTRAR_ALTLOC = True

# ---- Diagnóstico ----------------------------------------------------
# Desligado por padrão: o -help do Rosetta é uma tabela whose option names
# are split/wrapped across lines, so flags_suportadas() cannot recover the
# names reliably and reports EVERY flag as unknown. A check that always
# cries wolf is worse than no check. Set True only while changing the
# command line, and confirm any warning with:
#     ligand_dock -help | grep <flag>
VALIDAR_FLAGS = False
SEMENTE_BASE = 20260923    # base para os seeds; corridas distintas
                           # recebem seeds distintos e reprodutíveis

# =====================================================================
# --------------------- 2. CAMINHOS E DIRETÓRIOS ----------------------
# =====================================================================

ROSETTA_BIN_LIGAND   = "ligand_dock"                     # resolvido pelo PATH
ROSETTA_DB           = Path(os.environ["ROSETTA3_DB"])
NOME_LIGANTE         = "LG1"                             # código de 3 letras

WORK_DIR      = Path("/app/data")
PROTEINS_DIR  = WORK_DIR / "proteins"     # receptores .pdb (entrada)
LIGANDS_DIR   = WORK_DIR / "ligands"      # ligantes .sdf  (entrada)
PARAMS_DIR    = WORK_DIR / "params"       # .params e PDB do ligante
PREP_DIR      = WORK_DIR / "prep"         # receptores limpos + complexos iniciais
RESULTS_DIR   = WORK_DIR / "resultados"   # saída do Rosetta + ranking
LOGS_DIR      = WORK_DIR / "logs"         # stdout/stderr por corrida

for _d in (PROTEINS_DIR, LIGANDS_DIR, PARAMS_DIR, PREP_DIR, RESULTS_DIR, LOGS_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# Registro tag -> ponto de partida, preenchido em main() e consumido pelo
# diagnóstico de deslocamento no fim da execução.
CENTROS_PARTIDA = {}


# =====================================================================
# --------------------- 3. VALIDAÇÃO DE AMBIENTE ----------------------
# =====================================================================

def validar_ambiente():
    problemas = []

    if shutil.which(ROSETTA_BIN_LIGAND) is None:
        problemas.append(f"Binário não encontrado no PATH: {ROSETTA_BIN_LIGAND}")
    if not ROSETTA_DB.is_dir():
        problemas.append(f"Database do Rosetta não encontrado: {ROSETTA_DB}")
    if not PROTEINS_DIR.is_dir() or not list(PROTEINS_DIR.glob("*.pdb")):
        problemas.append(f"Nenhum .pdb em {PROTEINS_DIR}")
    if not LIGANDS_DIR.is_dir() or not list(LIGANDS_DIR.glob("*.sdf")):
        problemas.append(f"Nenhum .sdf em {LIGANDS_DIR}")
    if MODO_BUSCA not in ("GLOBAL", "LOCAL"):
        problemas.append(f"MODO_BUSCA inválido: {MODO_BUSCA}")
    if MODO_BUSCA == "LOCAL" and CENTRO_SITIO == (0.0, 0.0, 0.0):
        problemas.append("MODO_BUSCA='LOCAL' exige CENTRO_SITIO preenchido "
                         "com as coordenadas reais do bolsão.")

    if problemas:
        raise SystemExit("ERROS DE CONFIGURAÇÃO:\n  - " + "\n  - ".join(problemas))


def flags_suportadas():
    try:
        r = subprocess.run([str(ROSETTA_BIN_LIGAND), "-help"],
                           capture_output=True, text=True, timeout=120)
        texto = (r.stdout or "") + (r.stderr or "")
        # As flags aparecem no help como -namespace:sub:opcao
        return set(re.findall(r"(-[A-Za-z][\w:]*)", texto))
    except Exception as e:
        print(f"[aviso] não foi possível ler o -help do binário ({e}). "
              f"Seguindo sem validar flags.")
        return None


def checar_flags(cmd, suportadas):
    """Avisa sobre qualquer flag do comando ausente no -help."""
    if not suportadas:
        return
    desconhecidas = [t for t in cmd
                     if t.startswith("-") and not t[1:2].isdigit()
                     and t not in suportadas]
    if desconhecidas:
        print(f"[AVISO] flags não reconhecidas pela sua build: {desconhecidas}\n"
              f"        Confirme com: {ROSETTA_BIN_LIGAND} -help | grep <flag>")


# =====================================================================
# --------------------- 4. GEOMETRIA / PDB (utilitários) --------------
# =====================================================================
# Um PDB é formato de colunas fixas. Nunca faça split() em linha de PDB:
# campos podem colar entre si em estruturas grandes.
#   cols 13-16 nome do átomo | col 17 altloc | cols 18-20 resname
#   col 22 cadeia | cols 31-54 x,y,z | cols 77-78 elemento

def _ler_xyz(linha):
    """Extrai (x, y, z) de uma linha ATOM/HETATM por posição de coluna."""
    return (float(linha[30:38]), float(linha[38:46]), float(linha[46:54]))


def _escrever_xyz(linha, xyz):
    """O formato %8.3f é obrigatório: 8 caracteres, 3 decimais. Qualquer
    desvio quebra o parser do Rosetta.
    """
    x, y, z = xyz
    return f"{linha[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{linha[54:]}"


def centroide(pontos):
    """Centro geométrico (não ponderado por massa — suficiente aqui)."""
    n = len(pontos)
    return (sum(p[0] for p in pontos) / n,
            sum(p[1] for p in pontos) / n,
            sum(p[2] for p in pontos) / n)


def limpar_receptor(prot_path):
    """Prepara o receptor: filtra altlocs, trata HETATM, preserva cadeias."""
    saida = PREP_DIR / f"{prot_path.stem}_clean.pdb"
    cadeias_vistas = set()
    linhas = []

    with open(prot_path) as f:
        for linha in f:
            if linha.startswith("TER"):
                linhas.append("TER\n")
                continue
            if not linha.startswith(("ATOM", "HETATM")):
                continue

            # --- altloc: aceita branco ou o primeiro rótulo ("A") ---
            if FILTRAR_ALTLOC and linha[16] not in (" ", "A"):
                continue

            resname = linha[17:20].strip().upper()

            if linha.startswith("HETATM"):
                # Água é sempre removida: o Rosetta trata solvatação
                # implicitamente e águas explícitas atrapalham a grid.
                if resname in ("HOH", "WAT", "DOD"):
                    continue
                if resname not in HETATM_MANTER:
                    continue
                # Íon mantido precisa existir na database do Rosetta;
                # metais comuns existem, cofatores orgânicos exigiriam
                # um .params próprio.

            cadeia = linha[21]
            if cadeia == "X":
                # "X" fica reservada ao ligante. Remapeia para "W".
                linha = linha[:21] + "W" + linha[22:]
                cadeia = "W"
            cadeias_vistas.add(cadeia)
            linhas.append(linha)

    with open(saida, "w") as out:
        out.writelines(linhas)
        out.write("END\n")

    coords = [_ler_xyz(l) for l in linhas if l.startswith(("ATOM", "HETATM"))]
    print(f"  receptor limpo: {saida.name} "
          f"({len(coords)} átomos, cadeias {sorted(cadeias_vistas)})")

    if len(cadeias_vistas) > 1:
        print("  [nota] receptor multimérico — confirme se o sítio de "
              "interesse está numa interface antes de interpretar o score.")
    return saida, coords


# =====================================================================
# --------------------- 5. DETECÇÃO DE PONTOS DE PARTIDA --------------
# =====================================================================
#
#  POR QUE ESTA SEÇÃO FOI REESCRITA
#  --------------------------------
#  A versão anterior escolhia os pontos que tinham MAIS átomos num raio
#  de 8 A ("mais enterrados"). Isso está conceitualmente errado: o ponto
#  mais enterrado de uma proteína globular é o NÚCLEO HIDROFÓBICO, que
#  não é um sítio de ligação e sequer é acessível ao solvente.
#
#  Na prática, com N_PONTOS_PARTIDA = 5 e separação mínima de 10 A, os
#  cinco pontos "vencedores" cabiam todos dentro do mesmo núcleo. As
#  cinco corridas partiam essencialmente do mesmo lugar — daí o ligante
#  aparecer sempre na mesma posição mesmo no modo GLOBAL.
#
#  CRITÉRIO NOVO: um bolsão é um ponto CÔNCAVO MAS COM SAÍDA.
#  Medimos isso por lançamento de raios: de cada ponto candidato,
#  disparamos raios em direções distribuídas uniformemente sobre a
#  esfera e contamos quantos batem em átomos da proteína.
#
#      fração bloqueada ~ 1,00  -> núcleo interno (descartar)
#      fração bloqueada ~ 0,3   -> superfície plana/convexa (descartar)
#      fração bloqueada 0,6–0,95 -> concavidade com saída = BOLSÃO
#
#  É o mesmo princípio dos detectores de verdade (fpocket, POCASA,
#  LIGSITE). Continua sendo uma aproximação: se você tiver a saída de um
#  fpocket, use MODO_BUSCA="LOCAL" com aquele centro.
# =====================================================================

def _grade_espacial(coords, celula=8.0):
    """Indexa os átomos em células cúbicas para busca de vizinhos O(1).

    Sem isso, testar milhares de pontos candidatos contra milhares de
    átomos seria O(n*m). Com o lançamento de raios (30 direções x ~10
    passos por ponto) o custo sem indexação seria proibitivo.
    """
    grade = {}
    for c in coords:
        chave = (int(c[0] // celula), int(c[1] // celula), int(c[2] // celula))
        grade.setdefault(chave, []).append(c)
    return grade, celula


def _vizinhos(ponto, grade, celula):
    """Devolve os átomos das 27 células ao redor do ponto."""
    cx, cy, cz = (int(ponto[0] // celula),
                  int(ponto[1] // celula),
                  int(ponto[2] // celula))
    out = []
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for k in (-1, 0, 1):
                out.extend(grade.get((cx + i, cy + j, cz + k), ()))
    return out


def _direcoes_esfera(n=30):
    """Gera n direções quase uniformes sobre a esfera (espiral de Fibonacci).

    Distribuir bem as direções importa: se os raios forem só os 6 eixos
    cartesianos, a fração bloqueada vira muito ruidosa e bolsões rasos
    passam ou são rejeitados por acaso.
    """
    dirs = []
    phi = math.pi * (3.0 - math.sqrt(5.0))      # ângulo áureo
    for i in range(n):
        y = 1 - (i / float(n - 1)) * 2          # varre de +1 a -1
        raio = math.sqrt(max(0.0, 1 - y * y))
        theta = phi * i
        dirs.append((math.cos(theta) * raio, y, math.sin(theta) * raio))
    return dirs


DIRECOES = _direcoes_esfera(30)


def _fracao_bloqueada(p, grade, celula, alcance=14.0, passo=1.5, raio_atomo=2.7):
    """Fração dos raios saindo de `p` que colidem com a proteína.

    Caminha ao longo de cada direção em passos de `passo` angstroms até
    `alcance`. Se em algum passo houver átomo a menos de `raio_atomo`,
    o raio está bloqueado e paramos aquele raio.

    Interpretação do valor devolvido:
      1.0  = ponto completamente encerrado -> interior da proteína
      0.0  = ponto em campo aberto        -> solvente
      alto mas < 1 = concavidade com abertura -> bolsão
    """
    bloqueados = 0
    r2 = raio_atomo * raio_atomo
    for d in DIRECOES:
        t = passo
        while t <= alcance:
            q = (p[0] + d[0] * t, p[1] + d[1] * t, p[2] + d[2] * t)
            viz = _vizinhos(q, grade, celula)
            if any((q[0]-a[0])**2 + (q[1]-a[1])**2 + (q[2]-a[2])**2 <= r2
                   for a in viz):
                bloqueados += 1
                break
            t += passo
    return bloqueados / len(DIRECOES)


def detectar_cavidades(coords, n_pontos, etiqueta="receptor"):
    """Propõe centros de bolsão por concavidade acessível.

    Etapa 1 (barata): varre uma grade sobre a caixa da proteína e guarda
        pontos que estão FORA dos átomos mas em CONTATO com eles, com
        enterramento dentro de uma JANELA (não mais "quanto maior melhor").
    Etapa 2 (cara): lança raios só nos melhores candidatos da etapa 1 e
        mede a fração bloqueada, descartando núcleo e superfície plana.
    Etapa 3: desbaste guloso por separação mínima, para as N corridas não
        caírem todas na mesma cavidade.
    """
    PASSO = 1.5                 # resolução da varredura, em angstroms
    D_MIN, D_MAX = 3.3, 5.5     # distância ao átomo mais próximo
                                #  < 3,3 = dentro do volume de van der Waals
                                #  > 5,5 = longe demais para fazer contato
    RAIO_ENTERRAMENTO = 8.0
    # JANELA de enterramento. O teto é a correção central: sem ele, o
    # ranking premia o núcleo da proteína.
    ENTERRAMENTO_MIN, ENTERRAMENTO_MAX = 28, 70
    # Faixa de concavidade aceitável na etapa de raios.
    FRACAO_MIN, FRACAO_MAX = 0.58, 0.96
    SEPARACAO = 14.0            # separação mínima entre pontos aceitos.
                                # Era 10 A — pequeno demais, cabiam vários
                                # pontos dentro do mesmo bolsão/núcleo.
    MAX_CANDIDATOS_RAIOS = 400  # teto de custo da etapa 2

    grade, celula = _grade_espacial(coords, celula=RAIO_ENTERRAMENTO)
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]

    # ---------------- Etapa 1: triagem geométrica barata ----------------
    candidatos = []
    x = min(xs)
    while x <= max(xs):
        y = min(ys)
        while y <= max(ys):
            z = min(zs)
            while z <= max(zs):
                p = (x, y, z)
                viz = _vizinhos(p, grade, celula)
                if len(viz) >= ENTERRAMENTO_MIN:
                    d2_min = min((p[0]-a[0])**2 + (p[1]-a[1])**2 + (p[2]-a[2])**2
                                 for a in viz)
                    if D_MIN**2 <= d2_min <= D_MAX**2:
                        n_prox = sum(
                            1 for a in viz
                            if (p[0]-a[0])**2 + (p[1]-a[1])**2 + (p[2]-a[2])**2
                            <= RAIO_ENTERRAMENTO**2)
                        # JANELA, não máximo: descarta o núcleo já aqui
                        if ENTERRAMENTO_MIN <= n_prox <= ENTERRAMENTO_MAX:
                            candidatos.append((n_prox, p))
                z += PASSO
            y += PASSO
        x += PASSO

    if not candidatos:
        print(f"  [AVISO] nenhum ponto passou na triagem geométrica. "
              f"Usando o centroide — resultado não interpretável.")
        return [centroide(coords)]

    # ---------------- Etapa 2: concavidade por lançamento de raios -------
    candidatos.sort(key=lambda t: -t[0])
    bolsoes = []
    for n_prox, p in candidatos[:MAX_CANDIDATOS_RAIOS]:
        frac = _fracao_bloqueada(p, grade, celula)
        if FRACAO_MIN <= frac <= FRACAO_MAX:
            bolsoes.append((frac, n_prox, p))

    if not bolsoes:
        print("  [AVISO] nenhuma concavidade acessível encontrada. "
              "Possíveis causas: receptor muito pequeno, estrutura sem "
              "bolsão definido, ou parâmetros da janela apertados demais.")
        return [centroide(coords)]

    # Mais côncavo primeiro (dentro da faixa aceitável)
    bolsoes.sort(key=lambda t: -t[0])

    # ---------------- Etapa 3: desbaste por separação --------------------
    selecionados = []
    for frac, n_prox, p in bolsoes:
        if all((p[0]-q[0])**2 + (p[1]-q[1])**2 + (p[2]-q[2])**2 > SEPARACAO**2
               for _, _, q in selecionados):
            selecionados.append((frac, n_prox, p))
        if len(selecionados) >= n_pontos:
            break

    # ---------------- Relatório e arquivo de inspeção visual -------------
    print(f"  {len(selecionados)} bolsão(ões) candidato(s):")
    for i, (frac, n_prox, p) in enumerate(selecionados, start=1):
        print(f"    p{i:02d}  centro=({p[0]:7.2f},{p[1]:7.2f},{p[2]:7.2f})  "
              f"concavidade={frac:.2f}  vizinhos={n_prox}")

    # Distância entre os pontos — diagnóstico direto do problema que você
    # observou. Se todas as distâncias forem pequenas, os pontos estão na
    # mesma região e o "blind docking" é ilusório.
    if len(selecionados) > 1:
        dists = []
        for i in range(len(selecionados)):
            for j in range(i + 1, len(selecionados)):
                a, b = selecionados[i][2], selecionados[j][2]
                dists.append(math.dist(a, b))
        print(f"    separação entre pontos: min {min(dists):.1f} A | "
              f"max {max(dists):.1f} A")
        if min(dists) < SEPARACAO * 0.9:
            print("    [AVISO] pontos mais próximos que o esperado.")

    escrever_pontos_pdb(selecionados, etiqueta)
    return [p for _, _, p in selecionados]


def escrever_pontos_pdb(selecionados, etiqueta):
    """Grava os pontos de partida como pseudo-átomos, para inspeção visual.

    Este arquivo existe para você NÃO precisar confiar no script. Abra no
    PyMOL junto com o receptor limpo:

        load prep/<receptor>_clean.pdb
        load prep/<receptor>_pontos.pdb
        show spheres, resn PNT
        color red, resn PNT

    Se as esferas estiverem todas amontoadas, ou enfiadas no meio da
    proteína, o problema é a detecção de cavidade — não o docking.
    """
    caminho = PREP_DIR / f"{etiqueta}_pontos.pdb"
    with open(caminho, "w") as f:
        for i, (frac, n_prox, p) in enumerate(selecionados, start=1):
            # B-factor recebe a concavidade, para colorir por qualidade
            f.write(f"HETATM{i:5d}  O   PNT X{i:4d}    "
                    f"{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}"
                    f"  1.00{frac*100:6.2f}           O\n")
        f.write("END\n")
    print(f"    pontos gravados para inspeção: {caminho.name}")


# =====================================================================
# --------------------- 6. PREPARAÇÃO DO LIGANTE ----------------------
# =====================================================================

def preparar_ligante_params(lig_path):
    """Converte o SDF em .params + PDB compatível com Rosetta.

    SOBRE FLEXIBILIDADE DO LIGANTE — leia antes de confiar no resultado:
    a flag -docking:ligand:random_conformer só faz alguma coisa se o
    arquivo .params contiver uma linha PDB_ROTAMERS apontando para uma
    biblioteca de confôrmeros. Sem isso o docking é RÍGIDO e a flag é
    ignorada silenciosamente — o pior tipo de erro, porque não aparece
    em lugar nenhum do log.

    A função abaixo detecta a presença de PDB_ROTAMERS e avisa.
    """
    tag = lig_path.stem
    out_params = PARAMS_DIR / f"{tag}.params"
    out_pdb = PARAMS_DIR / f"{tag}_0001.pdb"

    # Cache: só reparametriza se faltar algum dos dois arquivos.
    # A versão anterior checava só o .params e devolvia um PDB que podia
    # não existir, causando falha obscura mais adiante.
    if not (out_params.exists() and out_pdb.exists()):
        print(f"  parametrizando {lig_path.name} ...")
        gerar_params(str(lig_path), NOME_LIGANTE, str(out_params))
        gerar_pdb(str(lig_path), NOME_LIGANTE, str(out_pdb), str(out_params))
        if not (out_params.exists() and out_pdb.exists()):
            raise RuntimeError(
                f"Parametrização falhou para {lig_path.name}.\n"
                f"Verifique se o SDF tem hidrogênios explícitos e coordenadas 3D.")

    # Diagnóstico de flexibilidade
    texto = out_params.read_text()
    tem_rotamers = "PDB_ROTAMERS" in texto
    if not tem_rotamers:
        print(f"  [AVISO] {out_params.name} sem PDB_ROTAMERS -> docking RÍGIDO. "
              f"Gere confôrmeros se o ligante tiver ligações rotacionáveis.")
    return out_params, out_pdb, tem_rotamers


def _matriz_rotacao_aleatoria(rng):
    """Matriz de rotação uniforme sobre SO(3), via quatérnio aleatório.

    POR QUE ISTO EXISTE: sem rotacionar, todas as corridas começam com o
    ligante na MESMA orientação herdada do SDF. Combinado com um
    protocolo que só refina localmente, isso produz poses quase
    idênticas entre pontos de partida — parte do sintoma de "o ligante
    ficou no mesmo lugar".

    O rng é semeado a partir do tag, então a rotação é ALEATÓRIA mas
    REPRODUTÍVEL: rodar de novo dá exatamente o mesmo complexo inicial.
    """
    u1, u2, u3 = rng.random(), rng.random(), rng.random()
    a = math.sqrt(1 - u1) * math.sin(2 * math.pi * u2)
    b = math.sqrt(1 - u1) * math.cos(2 * math.pi * u2)
    c = math.sqrt(u1) * math.sin(2 * math.pi * u3)
    d = math.sqrt(u1) * math.cos(2 * math.pi * u3)
    return (
        (1 - 2*(c*c + d*d), 2*(b*c - a*d),     2*(b*d + a*c)),
        (2*(b*c + a*d),     1 - 2*(b*b + d*d), 2*(c*d - a*b)),
        (2*(b*d - a*c),     2*(c*d + a*b),     1 - 2*(b*b + c*c)),
    )


def montar_complexo(receptor_pdb, lig_pdb, centro, tag, seed):
    """Monta o complexo inicial com o ligante TRANSLADADO para `centro`.

    Esta é a correção mais importante do script. Antes, o ligante era
    simplesmente concatenado nas coordenadas originais do SDF, que não
    têm relação nenhuma com o referencial do receptor. Aqui o centroide
    do ligante é movido para o ponto de partida escolhido, o que define
    onde a grid de pontuação do RosettaLigand será construída.

    A orientação inicial permanece a do SDF — a amostragem rotacional do
    protocolo se encarrega dela.
    """
    with open(lig_pdb) as f:
        linhas_lig = [l for l in f if l.startswith(("ATOM", "HETATM"))]

    coords_lig = [_ler_xyz(l) for l in linhas_lig]
    centro_lig = centroide(coords_lig)

    # --- rotação aleatória em torno do próprio centroide -----------------
    # Cada ponto de partida recebe uma orientação inicial diferente.
    R = _matriz_rotacao_aleatoria(random.Random(seed))
    coords_rot = []
    for (x, y, z) in coords_lig:
        vx, vy, vz = x - centro_lig[0], y - centro_lig[1], z - centro_lig[2]
        coords_rot.append((
            R[0][0]*vx + R[0][1]*vy + R[0][2]*vz,
            R[1][0]*vx + R[1][1]*vy + R[1][2]*vz,
            R[2][0]*vx + R[2][1]*vy + R[2][2]*vz,
        ))
    # Após a rotação o centroide está na origem, então a translação para o
    # ponto de partida é simplesmente somar `centro`.

    # NOME DO ARQUIVO — atenção, isto define o nome da saída do Rosetta.
    # O Rosetta nomeia cada pose como:  <out:prefix> + <nome do input> + _NNNN
    # Se o input se chamasse "{tag}_ini.pdb" E o prefixo fosse "{tag}_", o
    # resultado sairia como receptor_ligante_receptor_ligante_ini_0001.pdb.
    # Solução: o input carrega o tag puro e NÃO usamos -out:prefix.
    # Não há colisão porque entrada (prep/) e saída (results/) são pastas
    # diferentes.
    saida = PREP_DIR / f"{tag}.pdb"
    with open(saida, "w") as out:
        # --- receptor (cadeias originais preservadas) ---
        ultimo_res = 0
        with open(receptor_pdb) as f:
            for l in f:
                if l.startswith(("ATOM", "HETATM")):
                    out.write(l)
                    try:
                        ultimo_res = max(ultimo_res, int(l[22:26]))
                    except ValueError:
                        pass
        out.write("TER\n")
        # O ligante precisa de um resSeq que não colida com o da proteína.
        # O PDB do ligante traz "1", que se sobrepõe ao primeiro resíduo e
        # faz o Rosetta reportar todos os átomos do ligante como ausentes.
        res_lig = ultimo_res + 1
        # --- ligante, rotacionado e transladado, sempre na cadeia X ---
        for i, (l, (x, y, z)) in enumerate(zip(linhas_lig, coords_rot), start=1):
            nova = _escrever_xyz(l, (x + centro[0], y + centro[1], z + centro[2]))
            # Preserva colunas 7..21 (serial, nome, altloc, resname) e
            # substitui a cadeia (col 22) por "X". O resname ocupa 18..20,
            # então o corte em 21 é o que mantém o nome intacto.
            linha = "HETATM" + nova[6:21] + "X" + nova[22:]
            out.write(f"{linha[:22]}{res_lig:4d}{linha[26:]}")
        out.write("TER\nEND\n")
    return saida


# =====================================================================
# --------------------- 7. MONTAGEM DOS COMANDOS ----------------------
# =====================================================================

def montar_comando(complexo_pdb, params_file, tag, nstruct, seed, usar_conformeros):
    """Constrói a linha de comando do ligand_dock para uma corrida.

    Notas sobre as flags, já que a versão anterior tinha várias inventadas:
      -docking:ligand:protocol  precisa do namespace completo; "abbrev2" é
                                o protocolo padrão para ligantes pequenos
      -docking:uniform_trans    NÃO EXISTE. A opção real é uniform_trans_mag,
                                e ela pertence ao docking proteína-proteína.
                                No multi-start ela é desnecessária: a
                                diversidade vem dos pontos de partida.
      -run:jran                 seed. Corridas independentes precisam de
                                seeds distintos, senão você replica a mesma
                                trajetória e infla artificialmente a amostra.
    """
    cmd = [str(ROSETTA_BIN_LIGAND),
           "-database", str(ROSETTA_DB),
           "-s", str(complexo_pdb),
           "-extra_res_fa", str(params_file),   # obrigatório p/ molécula pequena
           "-nstruct", str(nstruct),
           "-docking:ligand:protocol", "abbrev2"]

    # --- diferença entre os modos ---
    if MODO_BUSCA == "GLOBAL":
        # Multi-start: a diversidade ENTRE corridas vem dos pontos de
        # partida. Mas dentro de cada corrida ainda é preciso amostrar
        # orientação, senão o ligante só relaxa a partir da pose inicial e
        # as poses saem quase idênticas — foi exatamente o que aconteceu na
        # versão anterior, que não tinha nenhuma flag de amostragem aqui.
        #
        # improve_orientation N: testa N orientações aleatórias no início e
        # fica com a melhor. É o passo de "encaixe grosseiro" do protocolo.
        cmd += ["-docking:ligand:improve_orientation", "100",
                "-docking:ligand:minimize_ligand", "true",
                "-docking:ligand:harmonic_torsions", "10"]
    else:
        # Refinamento local: perturbação pequena, mais minimização.
        cmd += ["-docking:ligand:minimize_ligand", "true",
                "-docking:ligand:minimize_backbone", "true",
                "-docking:ligand:harmonic_torsions", "10",
                "-docking:ligand:harmonic_Calphas", "0.2"]

    if usar_conformeros:
        cmd += ["-docking:ligand:random_conformer"]

    # --- tratamento da estrutura de entrada ---
    cmd += ["-run:preserve_header", "true",
            "-ignore_zero_occupancy", "false",
            "-ignore_unrecognized_res", "false",  # falhe alto se houver
                                                  # resíduo não reconhecido,
                                                  # em vez de sumir com ele
            "-ex1", "-ex2aro",                    # rotâmeros extras nas cadeias
                                                  # laterais da interface
            "-run:jran", str(seed)]

    # --- saída ---
    # Tudo cai em RESULTS_DIR. É o NOME DO INPUT que separa as corridas —
    # por isso as subpastas por complexo da versão anterior foram removidas
    # (eram criadas e nunca usadas, e faziam a coluna "complexo" do CSV
    # virar a constante "results").
    #
    # SEM -out:prefix de propósito: o input já se chama "{tag}.pdb", então
    # as poses saem como  receptor_ligante_p01_0001.pdb.
    # Acrescentar o prefixo aqui reintroduziria a duplicação do tag.
    cmd += ["-out:path:all", str(RESULTS_DIR),
            "-out:file:scorefile", f"{tag}_scores.sc",
            "-out:pdb", "true",
            "-overwrite"]
    return cmd


# =====================================================================
# --------------------- 8. EXECUÇÃO -----------------------------------
# =====================================================================

def executar(tarefa):
    """Roda uma corrida e devolve (tag, ok, mensagem).

    Usa ThreadPoolExecutor em vez de ProcessPoolExecutor: cada tarefa é
    apenas um subprocess bloqueante, ou seja, do ponto de vista do Python
    é I/O-bound. Processos aqui só adicionariam o custo de reimportar o
    módulo inteiro em cada worker (caro no macOS, que usa spawn).
    O paralelismo real é dos processos do Rosetta.
    """
    cmd, tag = tarefa
    log = LOGS_DIR / f"{tag}.log"
    try:
        with open(log, "w") as f:
            f.write("# " + " ".join(cmd) + "\n\n")
            f.flush()
            subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, check=True)
        return (tag, True, "ok")
    except subprocess.CalledProcessError as e:
        return (tag, False, f"código {e.returncode} — ver {log}")
    except KeyboardInterrupt:
        # Deixa a interrupção subir; o except nu da versão anterior
        # engolia Ctrl+C e a corrida virava imparável.
        raise
    except Exception as e:
        return (tag, False, f"{type(e).__name__}: {e}")


# =====================================================================
# --------------------- 9. ANÁLISE E RANKING --------------------------
# =====================================================================

def _chave_interface(cabecalho):
    """Localiza a coluna de energia de interface no score file.

    O app legado ligand_dock escreve 'interface_delta'. O
    InterfaceScoreCalculator do RosettaScripts escreve 'interface_delta_X',
    com o sufixo da cadeia do ligante. Aceitamos as duas formas para o
    script continuar válido se você migrar o protocolo depois.
    """
    for c in cabecalho:
        if c == "interface_delta" or c.startswith("interface_delta_"):
            return c
    return None


def _chave_rms(cabecalho):
    """RMSD do ligante, quando disponível."""
    for nome in ("ligand_rms_no_super_X", "ligand_rms_no_super", "rms"):
        if nome in cabecalho:
            return nome
    return None


def gerar_ranking():
    """Consolida todos os score files num CSV ordenado por afinidade.

    SOBRE A INTERPRETAÇÃO — o critério abaixo não é cosmético:

    Ordenar só por interface_delta é enganoso, porque poses com clash
    interno podem ter interface aparentemente favorável. O procedimento
    aceito na literatura do RosettaLigand é filtrar primeiro pelos
    melhores total_score (top ~5%) e só então ranquear por interface_delta
    dentro desse subconjunto.

    Valores POSITIVOS de interface_delta significam interface
    desfavorável — tipicamente clash ou desolvatação não compensada. Num
    blind docking multi-start é normal que pontos de partida ruins
    produzam poses positivas; elas são o sinal de que aquela cavidade não
    serve, não um erro do pipeline.
    """
    linhas = []
    for sc in sorted(RESULTS_DIR.glob("*_scores.sc")):
        cabecalho = None
        k_iface = k_rms = None
        with open(sc) as f:
            for linha in f:
                partes = linha.split()
                if not partes or partes[0] != "SCORE:":
                    continue
                if "total_score" in partes:
                    cabecalho = partes
                    k_iface = _chave_interface(cabecalho)
                    k_rms = _chave_rms(cabecalho)
                    continue
                if not cabecalho:
                    continue
                d = dict(zip(cabecalho, partes))

                # tag identifica a corrida; vem do NOME do arquivo, não da
                # pasta (a versão anterior usava score_file.parent.name, que
                # com -out:path:all era sempre "results")
                tag = sc.name.replace("_scores.sc", "")

                try:
                    total = float(d.get("total_score", "nan"))
                    iface = float(d.get(k_iface, "nan")) if k_iface else float("nan")
                except ValueError:
                    continue
                if math.isnan(total) or math.isnan(iface):
                    continue

                # tag = proteina_ligante_ponto (ver main); recupera o par
                partes_tag = tag.rsplit("_p", 1)
                par = partes_tag[0]
                ponto = partes_tag[1] if len(partes_tag) > 1 else "-"

                linhas.append({
                    "par": par,
                    "ponto_partida": ponto,
                    "pose_id": d.get("description", "?"),
                    "total_score": total,
                    "interface_delta": iface,
                    "rms_ligante": d.get(k_rms, "NA") if k_rms else "NA",
                })

    if not linhas:
        print("\nNenhum score file válido encontrado. Verifique os logs em "
              f"{LOGS_DIR}.")
        return

    # --- filtro por total_score, depois ranking por interface_delta ---
    linhas.sort(key=lambda r: r["total_score"])
    corte = max(1, int(len(linhas) * 0.05))
    elite = sorted(linhas[:corte], key=lambda r: r["interface_delta"])

    caminho_todos = RESULTS_DIR / "todas_as_poses.csv"
    caminho_rank = RESULTS_DIR / "ranking_ligantes.csv"
    colunas = ["par", "ponto_partida", "pose_id", "total_score",
               "interface_delta", "rms_ligante"]

    for caminho, dados in ((caminho_todos, linhas), (caminho_rank, elite)):
        with open(caminho, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=colunas)
            w.writeheader()
            w.writerows(dados)

    print(f"\n{len(linhas)} poses lidas | top {corte} por total_score ranqueadas")
    print(f"  todas as poses : {caminho_todos}")
    print(f"  ranking final  : {caminho_rank}")

    # --- melhor pose por par proteína-ligante ---
    print("\nMelhor pose por par (dentro do top 5% de total_score):")
    vistos = set()
    for r in elite:
        if r["par"] in vistos:
            continue
        vistos.add(r["par"])
        aviso = "  <-- interface desfavorável" if r["interface_delta"] > 0 else ""
        print(f"  {r['par']:<30} ponto {r['ponto_partida']:<3} "
              f"interface_delta={r['interface_delta']:8.2f} REU  "
              f"({r['pose_id']}){aviso}")


def verificar_deslocamento(centros):
    """Mede o quanto cada pose se afastou do seu ponto de partida.

    ESTE É O DIAGNÓSTICO QUE RESPONDE À PERGUNTA "o ligante ficou parado?".

    Para cada pose gerada, calcula o centroide do ligante (cadeia X) e
    compara com o centro de partida daquela corrida.

    Como ler o resultado:
      < 1 A   -> o ligante praticamente não se moveu. Ou o protocolo não
                 está amostrando, ou o ponto de partida já era um mínimo
                 muito profundo (raro).
      1–5 A   -> comportamento normal do RosettaLigand: refinamento local
                 dentro da grid. É o esperado.
      > 8 A   -> o ligante saiu da caixa original; verifique se não foi
                 empurrado para fora da proteína.

    E, o mais importante: se poses de PONTOS DIFERENTES convergirem para
    o mesmo centroide final, os pontos de partida não eram realmente
    distintos.
    """
    print("\n--- Diagnóstico de deslocamento (pose final vs ponto de partida) ---")
    resumo = {}
    for pdb in sorted(RESULTS_DIR.glob("*_[0-9][0-9][0-9][0-9].pdb")):
        # nome esperado: <tag>_0001.pdb
        tag = pdb.stem.rsplit("_", 1)[0]
        if tag not in centros:
            continue
        coords = []
        with open(pdb) as f:
            for l in f:
                if l.startswith(("ATOM", "HETATM")) and l[21] == "X":
                    coords.append(_ler_xyz(l))
        if not coords:
            continue
        c = centroide(coords)
        d = math.dist(c, centros[tag])
        resumo.setdefault(tag, []).append((d, c))

    if not resumo:
        print("  Nenhuma pose encontrada para medir.")
        return

    for tag in sorted(resumo):
        dados = resumo[tag]
        ds = [d for d, _ in dados]
        cfinal = centroide([c for _, c in dados])
        alerta = ""
        if max(ds) < 1.0:
            alerta = "   <-- PRATICAMENTE PARADO, verifique a amostragem"
        print(f"  {tag:<34} n={len(ds):<4} desloc {min(ds):5.1f}–{max(ds):5.1f} A"
              f"  centroide final ({cfinal[0]:7.2f},{cfinal[1]:7.2f},{cfinal[2]:7.2f})"
              f"{alerta}")

    # Convergência entre pontos distintos: sinal de que os pontos de
    # partida não eram independentes.
    finais = {t: centroide([c for _, c in v]) for t, v in resumo.items()}
    tags = sorted(finais)
    for i in range(len(tags)):
        for j in range(i + 1, len(tags)):
            if math.dist(finais[tags[i]], finais[tags[j]]) < 3.0:
                print(f"  [AVISO] {tags[i]} e {tags[j]} convergiram para o "
                      f"mesmo local — pontos de partida redundantes.")


# =====================================================================
# --------------------- 10. FLUXO PRINCIPAL ---------------------------
# =====================================================================

def main():
    validar_ambiente()
    suportadas = flags_suportadas() if VALIDAR_FLAGS else None

    proteinas = sorted(PROTEINS_DIR.glob("*.pdb"))
    ligantes = sorted(LIGANDS_DIR.glob("*.sdf"))
    fila = []

    print(f"\n=== MODO {MODO_BUSCA} | {len(proteinas)} receptor(es) x "
          f"{len(ligantes)} ligante(s) ===\n")

    # --- ligante no laço externo: parametriza uma vez, reaproveita ---
    for sdf in ligantes:
        try:
            params_file, lig_pdb, tem_conf = preparar_ligante_params(sdf)
        except Exception as e:
            print(f"[ERRO] ligante {sdf.name}: {e}")
            continue

        for prot in proteinas:
            print(f"\n[{prot.stem} + {sdf.stem}]")
            receptor, coords = limpar_receptor(prot)

            # --- define os pontos de partida ---
            if MODO_BUSCA == "GLOBAL":
                pontos = detectar_cavidades(coords, N_PONTOS_PARTIDA,
                                            etiqueta=prot.stem)
                nstruct = POSES_POR_PONTO
                if len(pontos) < N_PONTOS_PARTIDA:
                    print(f"  [nota] só {len(pontos)} bolsão(ões) passaram nos "
                          f"critérios (pedidos: {N_PONTOS_PARTIDA}). "
                          f"Afrouxe a janela em detectar_cavidades se precisar "
                          f"de mais pontos.")
            else:
                pontos = [CENTRO_SITIO]
                nstruct = POSES_LOCAL
                print(f"  sítio informado: {CENTRO_SITIO}")

            # --- uma corrida independente por ponto de partida ---
            for i, centro in enumerate(pontos, start=1):
                tag = f"{prot.stem}_{sdf.stem}_p{i:02d}"
                # Seed determinístico a partir do tag. hash() do Python é
                # randomizado a cada processo (PYTHONHASHSEED), o que tornava
                # as corridas NÃO reprodutíveis entre execuções — trocado por
                # uma soma estável dos bytes do tag.
                seed = SEMENTE_BASE + sum(ord(ch) * (k + 1)
                                          for k, ch in enumerate(tag)) % 100000
                complexo = montar_complexo(receptor, lig_pdb, centro, tag, seed)
                CENTROS_PARTIDA[tag] = centro
                cmd = montar_comando(complexo, params_file, tag,
                                     nstruct, seed, tem_conf)
                if suportadas:
                    checar_flags(cmd, suportadas)
                fila.append((cmd, tag))

    if not fila:
        raise SystemExit("Nada a executar.")

    total_poses = sum(int(c[c.index("-nstruct") + 1]) for c, _ in fila)
    print(f"\n=== {len(fila)} corridas | {total_poses} poses no total ===")
    print(f"Estimativa grosseira: {total_poses * 1.5 / MAX_NUCLEOS / 60:.1f} h "
          f"em {MAX_NUCLEOS} núcleos (~1,5 min/pose, sem contar throttling).")
    print("Ctrl+C interrompe; os score files já escritos permanecem válidos.\n")

    concluidas = 0
    with ThreadPoolExecutor(max_workers=MAX_NUCLEOS) as ex:
        futuros = [ex.submit(executar, t) for t in fila]
        for fut in as_completed(futuros):
            tag, ok, msg = fut.result()
            concluidas += 1
            marca = "ok " if ok else "FALHOU"
            print(f"  [{concluidas}/{len(fila)}] {tag}: {marca} {'' if ok else msg}")

    verificar_deslocamento(CENTROS_PARTIDA)
    gerar_ranking()

    # --- orientação de próximo passo, para o pipeline não parar aqui ---
    if MODO_BUSCA == "GLOBAL":
        print("\nPRÓXIMO PASSO: identifique o ponto de partida vencedor no "
              "ranking, pegue o centroide do ligante naquela pose e rode de "
              "novo com MODO_BUSCA='LOCAL' e POSES_LOCAL alto. A triagem "
              "cega localiza o sítio; ela não produz a pose final.")
    print("\nCONTROLE OBRIGATÓRIO antes de reportar qualquer resultado: "
          "redocking de um complexo cristalográfico do mesmo alvo. Se o "
          "pipeline não recupera a pose conhecida com RMSD < 2 A, os "
          "números acima não distinguem 'não liga' de 'pipeline quebrado'.")


if __name__ == "__main__":
    main()