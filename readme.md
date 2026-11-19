# Simulação BEM 2D de Gotas em Microcanais (Fortran)

Este diretório contém a implementação em Fortran 90 para simular o escoamento bidimensional de gotas newtonianas imiscíveis em regime de Stokes (`Re = 0`) através de microcanais planos e microcapilares constritos, utilizando o **Método dos Elementos de Contorno (BEM)** com elementos quadráticos contínuos e acompanhamento Lagrangeano da interface.

O código corresponde à implementação descrita no trabalho de graduação e no artigo presentes neste repositório:
- `Projeto_Graduacao.pdf` — *Aplicação do Método dos Elementos de Contorno ao Escoamento de Gotas em Canais Convergentes* (UnB, 2016)
- `1-s2.0-S0301932217304895-main.pdf` — *Flow of emulsion drops through a constricted microcapillary channel* (*International Journal of Multiphase Flow* 103, 2018, pp. 141–150)
- `Manual.txt` — guia original de preenchimento dos arquivos de entrada.

---

## 1. Arquitetura do Código Fortran

O fluxo de trabalho é dividido em dois programas independentes:

1. **`Create_Mesh_2D/` (`Mesh.f90` → `Mesh.exe`)**: pré-processador que lê a geometria do canal, as gotas e as condições de contorno de um arquivo `.txt`, gera a malha quadrática (3 nós por elemento) e grava o arquivo `output.txt` além de contornos iniciais `.plt`.
2. **`Processor/` (`Processing.f90` → `Gota.exe`)**: solver BEM transiente que lê `output.txt`, monta e resolve o sistema linear em cada passo de tempo, evolui a interface da gota e exporta séries temporais e perfis da gota em formato Tecplot `.plt`.

### 1.1. `Create_Mesh_2D/` (Gerador de malha)

| Arquivo | Função |
| --- | --- |
| `variables.f90` | Módulo `global`: armazena dimensões da malha (`npts`, `nseg`, `ndrop`, `nele`, `nnode`), coordenadas (`xnode`), conectividade (`conn`), condições de contorno (`bc`) e parâmetros físicos/temporais (`lambda`, `Ca`, `aH`, `dt`, `max_pos`, `max_it`). |
| `Mesh.f90` | Programa principal (`PROGRAM Mesh`). Executa em sequência `getfilename` → `readmesh` → `meshdisc` → `format_bc` → `output_file` → `tecplot_contour` e termina com `PAUSE`. |
| `getfilename.f90` | Lê o nome-base (sem extensão `.txt`) do caso especificado em `file_name.dat`. |
| `readmesh.f90` | Lê `<caso>.txt`, calcula `aH = 2.0*drop_geom(1,4)/(xpts(npts,3) - xpts(1,3))`, o número total de elementos `nele` e de nós `nnode = 2*nele`. |
| `alloc_memory1.f90` / `alloc_memory2.f90` | Alocam os vetores de entrada da geometria e os vetores finais da malha. |
| `meshdisc.f90`, `calc_maxdim.f90`, `max_dimension.f90` | Calculam a escala máxima do domínio (`max_dim`) e chamam a discretização nodal. |
| `create_nodes_conn.f90`, `calc_arc.f90`, `circle_center.f90` | Geram os nós e a conectividade quadrática (`conn`) ao longo de segmentos retos, arcos circulares e circunferências das gotas. Os nós de gota começam em `drop_first_last_node(:,2)`. |
| `format_bc.f90` | Converte as condições de contorno por segmento para os 3 nós de cada elemento: `0` (velocidade prescrita), `1` (tração prescrita), `2` (perfil parabólico normal desenvolvido `3/2 * V_m * (1 - (2d/L)^2)` em segmentos retos) e condição de interface livre/gota. |
| `output_file.f90` | Escreve o arquivo `output.txt` consumido pelo `Processor`. |
| `tecplot_contour.f90` | Gera `wall_contour.plt` e `drop_contour.plt` para inspeção visual da malha inicial. |
| `dad_1.txt`, `dad_2.txt`, `dad_3.txt` | Casos de exemplo: canal reto (`dad_1`), canal constrito `L=3` garganta `0.5` (`dad_2`) e canal constrito `L=6` garganta `0.25` (`dad_3`). |

### 1.2. `Processor/` (Solver BEM e evolução temporal)

O programa principal `Prossesing` em `Processing.f90` executa os seguintes passos:

1. **Inicialização**:
   - `variables.f90` (módulo `global`) e `variables_build_HGB.f90` (módulo `var_build_HGB`) definem o estado global e os buffers de integração numérica.
   - `read_input.f90` lê `output.txt` e chama `alloc_memory.f90` para alocar `xnode`, `conn`, `bc`, `H`, `G`, `B`, `matA`, `vector_b`, `sol`, `nd_vel` e `nd_trac`.
   - `first_drop_elem.f90` localiza `drop_f_elem`, o índice do primeiro elemento pertencente à interface da gota (`conn(el,5) == 1`).

2. **Loop temporal (`iterations = 1 .. max_it`, passo `dt`)**:
   Cada passo temporal executa duas subetapas (`imp_evol = 1, 2`) formando um esquema **preditor-corretor de 2ª ordem**:
   - **Montagem das matrizes (`building_HGB_sing.f90` — paralelizada com OpenMP)**:
     - Aloca buffers auxiliares (`alloc_memory_build_HGB.f90`) e gera pontos de quadratura de Gauss-Legendre regular (`npgr = 6`, `Gauss_Legendre.f90`) e Gauss logarítmica singular (`npgs = 10`, `Gauss_Logaritimic.f90`, `transform_eta_qsi.f90`).
     - Avalia funções de forma quadráticas (`shape_functions.f90`, `shape_functions_sing.f90`) e parâmetros geométricos como normal, tangente, jacobiano e curvatura nodal (`geometric_parameters.f90`, `calc_radius.f90`, `calc_curv_nodes.f90`, `calc_circle_curv.f90`, `curv_sign.f90`).
     - O laço externo sobre os elementos (`do el = 1, nele`) é paralelizado com OpenMP (`!$omp parallel do ... reduction(+:drop_area) schedule(static)`), mantendo buffers geométricos e matrizes locais privados por thread. Integra contribuições regulares em `regular_integral.f90` (com `mat_prod_2x2_2x6.f90`) e trata singularidades quando o nó fonte pertence ao elemento via separação singular + regular (`sing_integral.f90` + `sing_reg_integral.f90`).
     - Calcula os blocos diagonais da matriz `H` indiretamente pela soma negativa dos blocos fora da diagonal (`!$omp parallel do` sobre `m = 1, nnode`) e computa a área instantânea da gota (`calc_drop_area.f90`), guardando `drop_area_ini` na primeira avaliação.
   - **Aplicação das condições de contorno (`apply_bc2.f90`)**:
     - Combina `H`, `G` e `B` (que contém o termo capilar `~ 1/Ca`) formando o sistema linear `matA * sol = vector_b` (com auxílio de `prod_Axb.f90` paralelizado com OpenMP) e registra o mapeamento de incógnitas em `desl_pr`.
   - **Solução do sistema linear (`lapack_solver.f90` → `solve_lapack` / LAPACK `dgesv`)**:
     - Resolve `matA * sol = vector_b` por fatoração LU direta com pivoteamento parcial via **LAPACK (`dgesv`)**, acompanhada de `dgemm` paralelizado com OpenMP para acelerar a fatoração em múltiplos núcleos. (O antigo solver iterativo `Conj_grad.f90` permanece no repositório para referência.)
   - **Reorganização (`reorganize.f90`)**:
     - Separa o vetor solução e as condições prescritas em velocidades nodais (`nd_vel`) e trações nodais (`nd_trac`).
   - **Evolução da interface**:
     - `imp_evol = 1`: `evolution_euler.f90` salva `xnode_old` e `nd_vel_old` e faz o preditor explícito (`x = x_old + dt * v`).
     - `imp_evol = 2`: `evolution_euler_implicit.f90` aplica o corretor trapezoidal de 2ª ordem (`x = x_old + 0.5 * dt * (v + v_old)`).

3. **Diagnósticos, saída e realocação de nós**:
   - `flow_error.f90`: calcula o erro percentual de conservação de vazão (`flow_err`) integrando velocidades na entrada e na saída pela regra de Simpson.
   - `driver_pressure.f90`: calcula a diferença de tração normal média entre entrada e saída (`d_pressure`).
   - `tecplot_contour.f90`: grava `pressure.plt`, `flow_error.plt` e `area_error.plt` em todo passo, além de `wall_contour.plt` e `drop_t=XXXXX.plt` no passo 1 e a cada 10 passos.
   - `nodes_relocation.f90`: chamada duas vezes ao final de cada passo para reposicionar os nós intermediários (`nd = 2`) dos elementos da gota segundo Wrobel et al. (2009), evitando distorção excessiva da malha.
   - **Critério de parada**: encerra antes de `max_it` se `xnode(drop_first_last_node(1,2), 2) > max_pos`.

> **Arquivos legados / não utilizados no fluxo principal**: `apply_bc.f90`, `building_HGB.f90`, `gauss_elimination.f90`, `singular_integral.f90`, `create_HGB.f90`, `write_HGB.f90`, `write_matA.f90` e a subrotina `BCGS` em `Conj_grad.f90` permanecem no diretório por histórico/depuração, mas não são chamados por `Processing.f90`.

---

## 2. Como Compilar

Os passos abaixo foram validados no Windows com **GNU Fortran (`gfortran`) 16.2.0 (MSYS2 UCRT64)**.

### Cuidados importantes antes de compilar
- **Limpeza de artefatos antigos**: se existirem arquivos `.mod` ou `.o` antigos (por exemplo de compilações anteriores em Linux), remova-os antes de compilar:
  ```powershell
  Remove-Item -ErrorAction SilentlyContinue *.mod, *.o
  ```
- **Ordem dos módulos**: os arquivos que definem módulos (`variables.f90` e, no `Processor`, também `variables_build_HGB.f90`) precisam ser passados **primeiro** na linha de comando do `gfortran`. (O antigo `Processor\roda.txt` os listava no final, o que falha em uma compilação limpa.)
- **`PATH` no Windows (MSYS2)**: caso outras ferramentas no `PATH` interfiram na execução interna do `f951.exe`, coloque `C:\msys64\ucrt64\bin` no início do `PATH` na sessão do terminal antes de invocar o `gfortran`.

### 2.1. Compilar o gerador de malha (`Create_Mesh_2D`)

No PowerShell, dentro de `Create_Mesh_2D`:

```powershell
$env:PATH = "C:\msys64\ucrt64\bin;" + $env:PATH
cd ".\Create_Mesh_2D"
Remove-Item -ErrorAction SilentlyContinue *.mod, *.o

gfortran -O2 -static -ffree-line-length-512 `
  variables.f90 alloc_memory1.f90 alloc_memory2.f90 calc_arc.f90 `
  calc_maxdim.f90 circle_center.f90 create_nodes_conn.f90 format_bc.f90 `
  getfilename.f90 max_dimension.f90 meshdisc.f90 output_file.f90 `
  readmesh.f90 tecplot_contour.f90 Mesh.f90 -o Mesh.exe
```

### 2.2. Compilar o solver BEM (`Processor`)

No PowerShell, dentro de `Processor`:

```powershell
$env:PATH = "C:\msys64\ucrt64\bin;" + $env:PATH
cd ".\Processor"
Remove-Item -ErrorAction SilentlyContinue *.mod, *.o

gfortran -O3 -march=native -fopenmp -fdefault-real-8 -fdefault-double-8 -static -ffree-line-length-512 "-Wl,--stack,536870912" `
  variables.f90 variables_build_HGB.f90 alloc_memory.f90 alloc_memory_build_HGB.f90 `
  apply_bc.f90 apply_bc2.f90 building_HGB.f90 building_HGB_sing.f90 `
  calc_circle_curv.f90 calc_curv_nodes.f90 calc_drop_area.f90 calc_radius.f90 `
  Conj_grad.f90 lapack_solver.f90 curv_sign.f90 driver_pressure.f90 evolution_euler.f90 `
  evolution_euler_implicit.f90 first_drop_elem.f90 flow_error.f90 gauss_elimination.f90 `
  Gauss_Legendre.f90 Gauss_Logaritimic.f90 geometric_parameters.f90 mat_prod_2x2_2x6.f90 `
  nodes_relocation.f90 prod_Axb.f90 read_input.f90 regular_integral.f90 reorganize.f90 `
  shape_functions.f90 shape_functions_sing.f90 sing_integral.f90 sing_reg_integral.f90 `
  tecplot_contour.f90 transform_eta_qsi.f90 write_HGB.f90 write_matA.f90 Processing.f90 `
  -llapack -lblas -o Gota.exe
```

#### Significado das flags do `Processor`
- **`-fopenmp`**: habilita a paralelização multi-thread via OpenMP na montagem das matrizes `H`, `G` e `B` (`building_HGB_sing.f90`), no produto matriz-vetor (`prod_Axb.f90`) e na multiplicação de blocos `dgemm` usada pelo LAPACK.
- **`-llapack -lblas`**: vincula as bibliotecas LAPACK e BLAS para resolução direta do sistema linear denso via `dgesv` (`lapack_solver.f90`).
- **`-fdefault-real-8 -fdefault-double-8` (obrigatórias)**: promovem variáveis `real` para 64 bits (dupla precisão IEEE 754) mantendo `double precision` (`0.0d0`) também em 64 bits nativos (evitando promoção acidental para `REAL(16)` em software e garantindo compatibilidade binária com a ABI de 64 bits do LAPACK `dgesv`).
- **`"-Wl,--stack,536870912"` (obrigatória no Windows)**: aumenta a stack para 512 MB (entre aspas no PowerShell por causa das vírgulas). Arrays automáticos locais como `matB(2*nnode,6*nele)` em `apply_bc2.f90` excedem o limite padrão de 2 MB no Windows. *(No Linux, o padrão de stack costuma ser gerido via `ulimit -s unlimited` e `-static` é opcional.)*
- **`-static`**: embute as bibliotecas do compilador no `.exe`, permitindo rodar `Mesh.exe` e `Gota.exe` sem depender do `PATH` do MSYS2.
- **`-O3 -march=native`**: acelera significativamente a montagem das integrais e os produtos matriz-vetor (`matmul`/`dot_product`) no BiCG.
- **`-ffree-line-length-512`**: aceita linhas extensas presentes nos fontes Fortran 90.

> **Adaptações de portabilidade presentes neste repositório**:
> - `Create_Mesh_2D/readmesh.f90`: correção do rótulo `100 format(20(f8.5))` (espaço faltante em `100format`).
> - `Create_Mesh_2D/output_file.f90`: inclusão das vírgulas nos descritores de `format` nas linhas `101` e `118`.
> - `Processor/Conj_grad.f90`: remoção de `print*,r` dentro do loop iterativo e remoção da alocação de `matAt(N,N)` na stack em `Conj_grad`, utilizando `matmul` e `dot_product`.

---

## 3. Como Rodar

### Passo 1 — Configurar o caso de entrada em `Create_Mesh_2D/`

1. Crie ou escolha um arquivo `<nome>.txt` (por exemplo `dad_1.txt`, `dad_2.txt` ou `dad_3.txt`). A estrutura lida por `readmesh.f90` segue a ordem abaixo:
   1. Linha comentada / cabeçalho (ignorada) + linha com `npts nseg ndrop`;
   2. Linha comentada (ignorada) + linha com `dt max_pos max_it` *(passo de tempo, posição máxima `x` para parada da gota e número máximo de iterações)*;
   3. `npts` linhas com os vértices do contorno externo: `id x y`;
   4. `nseg` linhas definindo os segmentos: `id p_ini p_fim raio_ou_0 tipo` (`raio_ou_0 = 0` para reto, `> 0` centro à esquerda, `< 0` centro à direita; `tipo = 0` parede, `-1` entrada, `-2` saída);
   5. `nseg` linhas de condições de contorno: `id tipo_x val_x tipo_y val_y` (`0` velocidade, `1` tração, `2` perfil parabólico com velocidade média em `val_x`);
   6. `nseg` linhas com número de elementos por segmento: `id num_elem`;
   7. `drop_nele` (elementos por gota);
   8. `lambda` (razão de viscosidade $\mu_d / \mu$);
   9. `Ca` (número capilar);
   10. `ndrop` linhas de geometria das gotas: `id x_c y_c raio`.
2. Edite `Create_Mesh_2D/file_name.dat` contendo apenas o nome do arquivo **sem `.txt`** (ex.: `dad_2`).

### Passo 2 — Gerar a malha (`Mesh.exe`)

Dentro de `Create_Mesh_2D/`:

```powershell
cd ".\Create_Mesh_2D"
.\Mesh.exe
```

O programa termina com `PAUSE` (em terminal interativo basta pressionar `Enter`; para automatizar em PowerShell use `` "`n" | .\Mesh.exe ``). Serão gerados:
- `output.txt`: malha completa, conectividade, BCs e parâmetros para o solver;
- `wall_contour.plt` e `drop_contour.plt`: contornos iniciais da parede e da gota para conferência.

### Passo 3 — Transferir `output.txt` para `Processor/`

```powershell
Copy-Item ".\Create_Mesh_2D\output.txt" ".\Processor\output.txt"
```

> **Dica para testes rápidos**: antes de rodar simulações completas (por exemplo `5000`–`8000` iterações), você pode editar o terceiro número da 2ª linha numérica de `Processor\output.txt` (`max_it`, ex.: `50`) para validar o comportamento e o tempo por passo.

### Passo 4 — Executar o solver (`Gota.exe`)

Dentro de `Processor/`:

```powershell
cd ".\Processor"
.\Gota.exe
```

### Arquivos de saída gerados em `Processor/`

| Arquivo | Conteúdo | Frequência de atualização |
| --- | --- | --- |
| `pressure.plt` | Duas colunas: tempo `t` e queda de pressão média `d_pressure` entre entrada e saída. | Todo passo |
| `flow_error.plt` | Duas colunas: tempo `t` e erro relativo percentual da vazão (`(Q_out - Q_in)/Q_in * 100`). | Todo passo |
| `area_error.plt` | Duas colunas: tempo `t` e erro relativo percentual da área da gota (`(A - A_ini)/A_ini * 100`). | Todo passo |
| `wall_contour.plt` | Coordenadas `(x, y)` dos nós do canal. | Passo 1 e a cada 10 passos |
| `drop_t=<t_ms>.plt` | Coordenadas `(x, y)` do contorno da gota no instante `nint(t * 1000)` em milésimos adimensionais (ex.: `drop_t=  125.plt` para `t = 0.125`). | Passo 1 e a cada 10 passos |
