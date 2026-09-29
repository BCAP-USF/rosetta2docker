# Docking de ligantes com Rosetta

Este projeto roda o **Rosetta** dentro de um container Docker e calcula
**onde e com que força uma molécula pequena (ligante) se encaixa numa
proteína**.

É um serviço de *docking molecular*: você entrega a proteína e o ligante,
e o programa gera centenas de posições possíveis do ligante, ranqueia as
melhores e mostra qual é o encaixe mais provável.

---

## Índice

1. [O que você precisa ter instalado](#1-o-que-você-precisa-ter-instalado)
2. [Como colocar seus arquivos no lugar](#2-como-colocar-seus-arquivos-no-lugar)
3. [Como rodar](#3-como-rodar)
4. [Onde ficam os resultados](#4-onde-ficam-os-resultados)
5. [Como entender os resultados](#5-como-entender-os-resultados)
6. [Ajustando o programa](#6-ajustando-o-programa)
7. [Problemas comuns](#7-problemas-comuns)
8. [Para quem quer ir além](#8-para-quem-quer-ir-além)

---

## 1. O que você precisa ter instalado

**Só o Docker.** É o único programa necessário.

O Docker é o software que cria o "container" — uma caixinha isolada que já
vem com o Rosetta instalado dentro. Você não precisa instalar o Rosetta
nem Python na sua máquina.

Para verificar se o Docker está instalado, abra o terminal e digite:

```bash
docker --version
```

Se aparecer um número de versão, está tudo certo. Se der erro, baixe pelo
site <https://docs.docker.com/get-docker/>.

> **Observação:** a primeira vez que você rodar, o Docker vai baixar
> vários gigabytes de imagem. Isso é normal e acontece só uma vez. Pode
> levar alguns minutos.

---

## 2. Como colocar seus arquivos no lugar

Você precisa de **dois arquivos**:

| Arquivo | O que é | Onde colocar |
|---|---|---|
| `minha_proteina.pdb` | A estrutura da proteína | `data/proteins/` |
| `meu_ligante.sdf` | A molécula pequena que quer dockar | `data/ligands/` |

Seus arquivos ficam **na sua máquina**, dentro da pasta `data` deste
projeto. Para ver o caminho completo, abra o terminal **nesta pasta** e
digite `pwd`.

### Exemplo com a proteína RAR e o ligante Curcumina

```
data/
├── proteins/
│   └── RAR.pdb          ← estrutura da proteína
└── ligands/
    └── Curvumina.sdf    ← ligante
```

Se as pastas ainda não existirem (num clone novo do projeto, por exemplo),
crie-as uma vez:

```bash
mkdir -p data/proteins data/ligands
```

Depois é só arrastar seus arquivos para dentro delas.

> **Nomes de arquivo:** o programa ignora a diferença entre maiúsculas e
> minúsculas, mas não aceita espaços. Use nomes sem espaço, como
> `ligante_1.sdf`.

---

## 3. Como rodar

**Passo 1 — Ligue o container do Rosetta:**

```bash
docker compose up -d rosetta
```

Espere aparecer algo como `Started`. Rode este comando **sempre** antes de
rodar um cálculo, pois o container se desliga sozinho depois de um tempo
sem uso.

**Passo 2 — Rode o programa de docking:**

```bash
docker exec rosetta python3 /app/src/rosetta_lig_sdf.py
```

**Pronto.** Agora é só esperar.

O programa vai mostrar na tela o que está fazendo, algo como:

```
=== MODO GLOBAL | 1 receptor(es) x 1 ligante(s) ===

[RAR + Curvumina]
  receptor limpo: RAR_clean.pdb (3547 átomos, cadeias ['A'])
  5 bolsão(ões) candidato(s):
    p01  centro=( -14.12,  -4.62,  -9.59)  concavidade=0.93
    ...
```

> **Quanto tempo demora?** Depende do tamanho da proteína e da quantidade
> de posições pedidas. O exemplo deste projeto (uma proteína de ~3.500
> átomos, 25 posições) leva **poucos minutos** em um notebook comum. Para
> uma análise mais completa, veja a [seção 6](#6-ajustando-o-programa).

**Passo 3 — Desligue o container** (opcional, quando terminar):

```bash
docker compose stop rosetta
```

### Acompanhar o progresso de um cálculo longo

Se você fechou o terminal sem querer, o cálculo **continua rodando**.
Para ver como está, abra um novo terminal nesta pasta e digite:

```bash
docker exec rosetta tail -f /app/data/logs/rosetta_lig_sdf.log
```

Esse é o log geral (progresso e ranking final). Cada bolsão também grava
um log próprio, como `data/logs/RAR_Curvumina_p03.log`, com o detalhe
técnico do Rosetta.

Para sair dessa visualização sem parar o cálculo, pressione `Ctrl + C`.

### Parar um cálculo

```bash
docker compose stop rosetta
```

Nada se perde: o programa guarda o que já terminou e você pode rodar de
novo depois.

---

## 4. Onde ficam os resultados

Tudo é gravado na pasta `data/resultados/`, **na sua máquina**. Você pode
abrir os arquivos direto no seu sistema operacional, sem precisar de
nenhum programa extra.

```
data/resultados/
├── ranking_ligantes.csv      ← START HERE: a tabela com o resultado
├── todas_as_poses.csv        ← todas as posições, não só as melhores
├── RAR_Curvumina_p03_0001.pdb  ← a melhor posição, pronta para ver
└── RAR_Curvumina_p03_scores.sc  ← as notas de energia de cada posição
```

### O que é cada arquivo

| Arquivo | Para que serve |
|---|---|
| **`ranking_ligantes.csv`** | **Comece por aqui.** Lista os melhores encaixes, do melhor para o pior. É a resposta final. |
| `todas_as_poses.csv` | Todas as posições geradas, caso queira inspecionar alguma específica. |
| `*_0001.pdb` | Uma posição possível do ligante. É um arquivo 3D — você pode abrir num visualizador molecular (PyMOL, Chimera, Discovery Studio) para ver o encaixe. |
| `*_scores.sc` | As notas de energia calculadas pelo Rosetta. Você não precisa abrir. |

### O que o programa avisa no final

```
25 poses lidas | top 5 por total_score ranqueadas
  todas as poses : /app/data/resultados/todas_as_poses.csv
  ranking final  : /app/data/resultados/ranking_ligantes.csv

Melhor pose por par (dentro do top 5% de total_score):
  RAR_Curvumina   ponto 03  interface_delta=   -3.23 REU  (RAR_Curvumina_p03_0001)
```

O nome entre parênteses é o **arquivo da melhor posição**. É só abrir
esse arquivo `.pdb` para ver.

---

## 5. Como entender os resultados

Abra o arquivo `ranking_ligantes.csv`. Ele tem estas colunas:

| Coluna | O que significa |
|---|---|
| `par` | Qual combinação proteína + ligante é esta linha |
| `ponto_partida` | De qual bolsão da proteína o cálculo começou |
| `pose_id` | **O nome do arquivo `.pdb` com esta posição** |
| `total_score` | Nota geral da posição. **Quanto menor, melhor.** |
| `interface_delta` | Quão bem o ligante se encaixa. **Quanto mais negativo, melhor.** |

### Regras rápidas para interpretar

- **`interface_delta` negativo** = bom encaixe. Os números mais negativos
  (ex.: `-3.23`) são os melhores.
- **`interface_delta` positivo** = encaixe ruim. Aquele bolsão provavelmente
  não serve para este ligante. Isso é um resultado válido, não um erro.
- **Comece sempre pela primeira linha.** Ela já é a melhor posição.

### ⚠️ Leia antes de confiar no resultado

> **Este programa encontra *onde* o ligante prefere se encaixar. Ele não
> garante que o encaixe esteja certo.**

Antes de usar qualquer resultado em pesquisa ou publicação, faça o
**controle de redocking**: rode o mesmo pipeline com um ligante que já
tenha posição experimental conhecida na mesma proteína. Se o programa
recuperar a posição conhecida com erro menor que 2 Å (RMSD), o pipeline
está confiável. Se não recuperar, os resultados não distinguem "não liga"
de "programa com defeito".

Sem esse controle, um número bonito não significa nada.

---

## 6. Ajustando o programa

Todas as opções ficam no **começo** do arquivo
`src/rosetta_lig_sdf.py`, numa seção chamada "PAINEL DE CONTROLE". É o
único lugar que você precisa editar.

### Quantas posições gerar

```python
N_PONTOS_PARTIDA = 5       # quantos bolsões testar
POSES_POR_PONTO  = 5       # posições por bolsão
```

O total é o **produto** desses dois números. O exemplo acima gera
5 × 5 = **25 posições** e leva poucos minutos.

Para uma busca mais completa, aumente os dois (por exemplo, 12 × 25 = 300
posições, algumas horas).

### Refinar a partir de um bolsão específico

A busca cega (acima) serve para **descobrir** o bolsão. Depois de
descobrir, o melhor é gerar muitas posições **só naquele bolsão**:

```python
MODO_BUSCA = "LOCAL"
CENTRO_SITIO = (23.375, -6.125, -5.094)   # o ponto vencedor
POSES_LOCAL = 300
```

Para achar as coordenadas do ponto vencedor, procure no log a linha do
bolsão com o mesmo número:

```
    p03  centro=(  23.38,  -6.12,  -5.09)  concavidade=0.90  vizinhos=51
```

Copie os três números para `CENTRO_SITIO` (o `ranking_ligantes.csv` mostra
apenas `03`, não as coordenadas):

```bash
docker exec rosetta grep centro /app/data/logs/rosetta_lig_sdf.log
```

### Cuidade do seu computador

```python
MAX_NUCLEOS = 4
```

Se o seu notebook ficar muito lento ou aquecer, diminua para `3`. Notebooks
sem ventoinha (MacBook Air M1, por exemplo) sofrem *throttling*
térmico quando usam todos os núcleos por muitas horas.

### Ions metálicos importantes

```python
HETATM_MANTER = {"ZN", "MG", "CA", "MN", "FE", "CU", "NA", "K"}
```

O programa descarta da proteína tudo que não for proteína — incluindo
águas e íons. Se o seu sítio tem um **zinco catalítico** e ele não estiver
na lista, o resultado será **fisicamente sem sentido**. Confira os íons da
sua estrutura e adicione o que faltar.

---

## 7. Problemas comuns

### `Nenhum .pdb em .../proteins` ou `Nenhum .sdf em .../ligands`

O programa não achou seus arquivos. Confira se estão dentro de
`data/proteins/` e `data/ligands/` — e não soltos na pasta `data/`.

### `Cannot connect to the Docker daemon` / `docker exec` falha

O container não está rodando. Ligue antes:

```bash
docker compose up -d rosetta
```

### `[AVISO] flags não reconhecidas pela sua build`

**Pode ignorar.** Esse aviso é um falso positivo conhecido: o
`-help` do Rosetta é uma tabela que quebra os nomes das opções em várias
linhas, e o verificador não consegue lê-los. Todas as opções usadas são
válidas e o cálculo roda normalmente.

Se quiser confirmar uma opção específica:

```bash
docker exec rosetta ligand_dock -help | grep extra_res_fa
```

### `[AVISO] docking RÍGIDO` / ligante sem conformeros

O ligante foi tratado como uma peça rígida, sem liberdade para girar. Se a
molécula tem partes flexíveis (ponteiras, cadeias laterais), os
resultados representam apenas **uma** das possíveis formas dela.

### `[ERRO] Unrecognized residue`

Bug já corrigido nesta versão, mas pode reaparecer se você editar a
conversão SDF. O `.pdb` do ligante e o arquivo `.params` precisam ter
exatamente o mesmo nome de resíduo (3 letras).

### `too many tries in fill_missing_atoms`

O ligante e o `.params` não estão batendo. Na prática isso quase sempre
significa que o resíduo do ligante se sobrepõe ao da proteína no arquivo
`.pdb` gerado. Apague `data/params/` e `data/prep/` e rode de novo.

---

## 8. Para quem quer ir além

### Estrutura do projeto

```
.
├── docker-compose.yml     # configuração do container
├── src/                   # o código
│   ├── rosetta_lig_sdf.py   # programa principal de docking
│   └── molfile_to_params.py  # converte o ligante para formato do Rosetta
└── data/                  # SEUS ARQUIVOS E RESULTADOS
    ├── proteins/          # ← coloque aqui as proteínas (.pdb)
    ├── ligands/           # ← coloque aqui os ligantes (.sdf)
    ├── params/            # (gerado) arquivos do ligante convertidos
    ├── prep/              # (gerado) complexos prontos para o Rosetta
    ├── resultados/        # (gerado) as poses e tabelas
    └── logs/              # (gerado) o registro de cada execução
```

Note que os dados ficam **fora** do container. É por isso que você não
perde nada ao desligá-lo, e pode usar os resultados com qualquer outro
programa.

### Outros serviços disponíveis

O `docker-compose.yml` define outros containers, se precisar:

| Serviço | Para que serve |
|---|---|
| `rosetta` | O que você usa neste guia |
| `rosetta-mpi` | Para cálculos paralelos com `mpirun` |
| `rosetta-jupyter` | Notebooks Python (roda em <http://localhost:8888>) |
| `rosetta-ml` | Machine learning (PyRosetta, TensorFlow) |

### Sobre o Rosetta

O Rosetta é um software livre de simulação molecular da Universidade de
Washington, usado amplamente em estrutura de proteínas e desenho de fármacos.
Documentação: <https://rosettacommons.org>

> **Uso comercial:** a licença do PyRosetta (a versão Python) **exige
> compra** para uso comercial. Ver <https://license.rosettacommons.org>.
