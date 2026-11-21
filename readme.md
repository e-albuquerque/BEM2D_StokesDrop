# Simulação BEM 2D de Gotas em Microcanais (Gmsh + Python + Fortran)

Este diretório contém a implementação híbrida **Python (Gmsh + `meshio`) + Fortran 90 (OpenMP + LAPACK)** para simular o escoamento bidimensional de gotas newtonianas imiscíveis em regime de Stokes (`Re = 0`) através de microcanais planos e microcapilares constritos, utilizando o **Método dos Elementos de Contorno (BEM)** com elementos quadráticos contínuos de 3 nós e acompanhamento Lagrangeano da interface.

O código corresponde à implementação descrita no trabalho de graduação e no artigo presentes neste repositório:
- `Projeto_Graduacao.pdf` — *Aplicação do Método dos Elementos de Contorno ao Escoamento de Gotas em Canais Convergentes* (UnB, 2016)
- `1-s2.0-S0301932217304895-main.pdf` — *Flow of emulsion drops through a constricted microcapillary channel* (*International Journal of Multiphase Flow* 103, 2018, pp. 141–150)

---

## 1. Arquitetura do Projeto

O fluxo de trabalho segue o mesmo padrão modular do projeto `BEM2D_elastostatic` e é dividido em duas etapas principais:

1. **Pré-processador Python + Gmsh (`input_data.py` e `preproc.py`)**:
   - Gera a malha quadrática de 2ª ordem (`line3` no contorno e `triangle6` no domínio) a partir dos arquivos de geometria do Gmsh (`dad_1.geo`, `dad_2.geo`, `dad_3.geo`);
   - Lê o arquivo `.msh` com `meshio`, encadeia os elementos de contorno do canal (no sentido anti-horário) e da interface da gota (no sentido horário a partir de $\theta = 0$), aplica as condições de contorno (incluindo o perfil parabólico de Poiseuille na entrada), plota a malha e as condições de contorno (`geometry.png`) e grava o arquivo `Processor/output.txt`.
2. **Solver BEM Transiente em Fortran 90 (`Processor/Processing.f90` → `Processor/Gota.exe`)**:
   - Lê `Processor/output.txt`, monta as matrizes de influência `H`, `G` e `B` em paralelo com **OpenMP**, resolve o sistema linear denso via **LAPACK (`dgesv`)** em cada subetapa do esquema preditor-corretor de 2ª ordem, evolui a interface da gota e exporta séries temporais e perfis em arquivos de dados `.dat`.
3. **Pós-processador e Animação em Python (`posproc.py` e `animate_drop.py`)**:
   - `posproc.py` lê os arquivos `.dat` exportados pelo solver Fortran (`pressure.dat`, `flow_error.dat`, `area_error.dat`, `wall_contour.dat`, `drop_t=*.dat`), imprime o resumo de validação numérica e gera os gráficos com `matplotlib` na pasta `examples/`;
   - `animate_drop.py` gera um vídeo de animação (`examples/drop_animation.mp4`, com opção ou fallback para `.gif`) mostrando a gota atravessando o microcapilar sincronizada com a curva de pressão motriz $\Delta P(t)$.

### 1.1. Pré-processamento, Pós-processamento e Animação em Python (`input_data.py`, `preproc.py`, `posproc.py`, `animate_drop.py` e `*.geo`)

| Arquivo | Função |
| --- | --- |
| `dad_1.geo` | Geometria Gmsh do Problema 1: canal plano reto ($L = 10H$, $H = 1$) com 1 gota ($R = 0,25$, $a/H = 0,50$), `470` elementos quadráticos (`940` nós). |
| `dad_2.geo` | Geometria Gmsh do Problema 2: microcapilar constrito ($L = 3H$, razão de contração $H/h = 2$) com 1 gota ($R = 0,30$, $a/H = 0,60$), `220` elementos quadráticos (`440` nós). |
| `dad_3.geo` | Geometria Gmsh do Problema 3: microcapilar longo com constrição severa ($L = 6H$, $H/h = 4$) com 1 gota ($R = 0,15$, $a/H = 0,30$), `820` elementos quadráticos (`1640` nós). |
| `input_data.py` | Define a função `input_data(problem=2)`, retornando os parâmetros físicos ($\lambda$, $\text{Ca}$, $a/H$), parâmetros temporais ($\Delta t$, `max_pos`, `max_it`), o nome do arquivo `.geo`/`.msh` e o dicionário `bound_cond` associando as `Physical Curve` do Gmsh (`"wall"`, `"outlet"`, `"inlet"`, `"drop_1"`) às condições de contorno. |
| `preproc.py` | Executa o Gmsh (`ensure_mesh_exists`), lê `<caso>.msh` via `meshio` (`compute_inodes`), calcula as condições de contorno nodais (`mount_bcs`), gera o gráfico da malha e condições de contorno (`show_geometry` → `geometry.png`) e escreve `Processor/output.txt` (`create_input_file`). |
| `posproc.py` | Lê `output.txt` (`read_solver_input`) e os arquivos `.dat` do solver (`read_solver_outputs`), exibe o resumo quantitativo de validação (`print_validation_summary`) e plota as figuras de evolução da gota, pressão motriz, erros de conservação e painel resumo 2x2 (`run_postprocessor`). |
| `animate_drop.py` | Gera animação em vídeo MP4 (via FFmpeg / `imageio-ffmpeg`) ou GIF (via Pillow) da gota se deformando ao longo do canal sincronizada com a evolução da queda de pressão $\Delta P(t)$ (`examples/drop_animation.mp4`). |

#### Tipos de condições de contorno em `input_data.py` (`bound_cond`)
- `type_x, type_y = 0`: velocidade prescrita nas direções globais $x$ ou $y$ ($u_x, u_y$);
- `type_x, type_y = 1`: tração superficial prescrita nas direções globais $x$ ou $y$ ($t_x, t_y$);
- `type_x = 2`: **perfil parabólico de Poiseuille desenvolvido** normal ao segmento reto de entrada, com velocidade média $V_m$ informada em `value_x`:
  $$u_n(d) = \frac{3}{2} V_m \left[1 - \left(\frac{2d}{H}\right)^2\right]$$
  decomposto automaticamente em componentes globais $(u_x, u_y)$ com `type = 0` nos 3 nós de cada elemento quadrático;
- `type_x, type_y = 3`: interface deformável da gota ($\Gamma_d$, `elem_type = 1`).

### 1.2. `Processor/` (Solver BEM em Fortran 90 e evolução temporal)

O programa principal `Prossesing` em `Processing.f90` executa os seguintes passos:

1. **Inicialização**:
   - `variables.f90` (módulo `global`) e `variables_build_HGB.f90` (módulo `var_build_HGB`) definem o estado global e os buffers de integração numérica.
   - `read_input.f90` lê `output.txt` (gerado por `preproc.py`) e chama `alloc_memory.f90` para alocar `xnode`, `conn`, `bc`, `H`, `G`, `B`, `matA`, `vector_b`, `sol`, `nd_vel` e `nd_trac`.
   - `first_drop_elem.f90` localiza `drop_f_elem`, o índice do primeiro elemento pertencente à interface da gota (`conn(el,5) == 1`).

2. **Loop temporal (`iterations = 1 .. max_it`, passo `dt`)**:
   Cada passo temporal executa duas subetapas (`imp_evol = 1, 2`) formando um esquema **preditor-corretor de 2ª ordem**:
   - **Montagem das matrizes (`building_HGB_sing.f90` — paralelizada com OpenMP e cache parede-parede)**:
     - Aloca buffers auxiliares (`alloc_memory_build_HGB.f90`) e calcula uma única vez no passo inicial os pontos de quadratura de Gauss-Legendre regular (`npgr = 6`, `Gauss_Legendre.f90`), Gauss logarítmica singular (`npgs = 10`, `Gauss_Logaritimic.f90`, `transform_eta_qsi.f90`) e as funções de forma quadráticas (`shape_functions.f90`, `shape_functions_sing.f90`).
     - **Cache dos blocos estacionários parede-parede (`H_wall_cache` e `G_wall_cache`)**: como os nós e elementos das paredes do canal, entrada e saída são fixos no tempo, suas integrais mútuas (que correspondem a 67% dos pares elemento-nó em `dad_2` e 90,5% em `dad_3`) são calculadas apenas na primeira montagem e restauradas em memória nos passos seguintes, integrando nos elementos de parede apenas os nós fontes pertencentes à gota em movimento.
     - O laço externo sobre os elementos (`do el = 1, nele`) é paralelizado com OpenMP (`!$omp parallel do ... reduction(+:drop_area) schedule(dynamic, 4)`), mantendo buffers geométricos e matrizes locais privados por thread. Integra contribuições regulares de forma direta por blocos $2\times 2$ em `regular_integral.f90` e trata singularidades via separação singular + regular (`sing_integral.f90` + `sing_reg_integral.f90`).
     - Calcula os blocos diagonais da matriz `H` indiretamente pela soma negativa dos blocos fora da diagonal (`!$omp parallel do` sobre `m = 1, nnode`) e computa a área instantânea da gota (`calc_drop_area.f90`), guardando `drop_area_ini` na primeira avaliação.
   - **Aplicação das condições de contorno (`apply_bc2.f90` — montagem direta sem cópia de `G`)**:
     - Constrói o mapeamento topológico `desl_pr` uma única vez no primeiro passo e monta `matA` e `vector_b` diretamente a partir de `H`, `G` e `B` sem alocar ou copiar matrizes auxiliares de tamanho $2N \times 6N_e$, multiplicando apenas as colunas com condições de contorno prescritas não-nulas (`bc_value /= 0`).
   - **Solução do sistema linear (`lapack_solver.f90` → `solve_lapack` / LAPACK `dgesv`)**:
     - Resolve `matA * sol = vector_b` por fatoração LU direta com pivoteamento parcial via **LAPACK (`dgesv`)**, acompanhada de `dgemm` paralelizado com OpenMP para acelerar a fatoração em múltiplos núcleos.
   - **Reorganização (`reorganize.f90`)**:
     - Separa o vetor solução e as condições prescritas em velocidades nodais (`nd_vel`) e trações nodais (`nd_trac`).
   - **Evolução da interface**:
     - `imp_evol = 1`: `evolution_euler.f90` salva `xnode_old` e `nd_vel_old` e faz o preditor explícito ($\mathbf{x}^* = \mathbf{x}^n + \Delta t\,\mathbf{u}^n$).
     - `imp_evol = 2`: `evolution_euler_implicit.f90` aplica o corretor trapezoidal de 2ª ordem ($\mathbf{x}^{n+1} = \mathbf{x}^n + \frac{1}{2}\Delta t\,(\mathbf{u}^n + \mathbf{u}^*)$).

3. **Diagnósticos, saída e realocação de nós**:
   - `flow_error.f90`: calcula o erro percentual de conservação de vazão (`flow_err`) integrando velocidades na entrada e na saída pela regra de Simpson.
   - `driver_pressure.f90`: calcula a diferença de tração normal média entre entrada e saída (`d_pressure`).
   - `write_output.f90`: grava `pressure.dat`, `flow_error.dat` e `area_error.dat` em todo passo, além de `wall_contour.dat` e `drop_t=XXXXX.dat` no passo 1 e a cada 10 passos.
   - `nodes_relocation.f90`: chamada duas vezes ao final de cada passo para reposicionar os nós dos elementos da gota segundo Wrobel et al. (2009), evitando distorção excessiva da malha.
   - **Critério de parada**: encerra antes de `max_it` se `xnode(drop_first_last_node(1,2), 2) > max_pos`.

---

## 2. Como Compilar o Solver Fortran (`Processor/Gota.exe` / `Processor/Gota`)

O projeto inclui um **`Makefile`** na raiz que detecta automaticamente o sistema operacional:
- **No Windows (`OS == Windows_NT`)**: gera `Processor/Gota.exe` com as flags `-static -Wl,--stack,536870912`;
- **No Linux / macOS**: gera `Processor/Gota` vinculando dinamicamente `-fopenmp -llapack -lblas`.
- Todos os objetos (`.o`) e módulos (`.mod`) são armazenados separadamente em `build/`.

### 2.1. Compilação via `Makefile` (Recomendada)

- **No Windows (PowerShell com MSYS2 UCRT64)**:
  ```powershell
  $env:PATH = "C:\msys64\ucrt64\bin;" + $env:PATH
  mingw32-make clean
  mingw32-make -j4
  ```
- **No Linux / Google Colab / macOS**:
  ```bash
  make clean
  make -j4
  ```

### 2.2. Compilação manual direta com `gfortran` (Opcional)

Caso prefira compilar manualmente sem `make`, execute dentro de `Processor/`:

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
  write_output.f90 transform_eta_qsi.f90 write_HGB.f90 write_matA.f90 Processing.f90 `
  -llapack -lblas -o Gota.exe
```

### Significado das flags do `Processor`
- **`-fopenmp`**: habilita a paralelização multi-thread via OpenMP na montagem das matrizes `H`, `G` e `B` (`building_HGB_sing.f90`), no produto matriz-vetor (`prod_Axb.f90`) e na multiplicação de blocos `dgemm` usada pelo LAPACK.
- **`-llapack -lblas`**: vincula as bibliotecas LAPACK e BLAS para resolução direta do sistema linear denso via `dgesv` (`lapack_solver.f90`).
- **`-fdefault-real-8 -fdefault-double-8` (obrigatórias)**: promovem variáveis `real` para 64 bits (dupla precisão IEEE 754) mantendo `double precision` (`0.0d0`) também em 64 bits nativos, compatíveis com a ABI do LAPACK `dgesv`.
- **`"-Wl,--stack,536870912"` (obrigatória no Windows)**: aumenta a stack para 512 MB para acomodar arrays automáticos locais como `matB(2*nnode,6*nele)` em `apply_bc2.f90`.
- **`-static`**: embute as bibliotecas do compilador no `.exe`, permitindo rodar `Gota.exe` sem depender do `PATH` do MSYS2.

---

## 3. Como Rodar

### Passo 1 — Instalar dependências Python (se necessário)

```powershell
pip install numpy matplotlib meshio gmsh imageio-ffmpeg
```

### Passo 2 — Executar o pré-processador Python (`preproc.py`)

Na raiz do projeto, execute `preproc.py` escolhendo o problema (`-p 1` para `dad_1`, `-p 2` para `dad_2` [padrão] ou `-p 3` para `dad_3`):

```powershell
python preproc.py -p 2 -o Processor/output.txt --plot geometry.png
```

O script `preproc.py` irá:
1. Gerar automaticamente o arquivo de malha quadrática `dad_2.msh` a partir de `dad_2.geo` usando o Gmsh;
2. Ler `dad_2.msh` via `meshio`, aplicar as condições de contorno definidas em `input_data.py` e salvar a figura da malha e condições de contorno em `geometry.png`;
3. Gravar o arquivo `Processor/output.txt` pronto para o solver Fortran.

### Passo 3 — Executar o solver Fortran (`Processor/Gota.exe`)

Dentro de `Processor/`:

```powershell
cd ".\Processor"
.\Gota.exe
```

### Passo 4 — Executar o pós-processador (`posproc.py`) e gerar animação em vídeo (`animate_drop.py`)

Na raiz do projeto, execute `posproc.py` para imprimir o resumo de validação e gerar todas as figuras estáticas em `examples/`:

```powershell
python posproc.py
```

Para gerar o vídeo de animação (`examples/drop_animation.mp4`) da gota atravessando o canal sincronizada com $\Delta P(t)$:

```powershell
python animate_drop.py
```

Opções disponíveis em `animate_drop.py`:
- `--fps <int>`: quadros por segundo (padrão: `15`);
- `--dpi <int>`: resolução em DPI (padrão: `150`);
- `--gif`: exporta também uma versão em GIF além do MP4;
- `-r <dir>` / `-o <arquivo>`: permite escolher outra pasta de resultados ou nome de saída.

### Arquivos de saída gerados em `Processor/`

| Arquivo | Conteúdo | Frequência de atualização |
| --- | --- | --- |
| `pressure.dat` | Duas colunas: tempo `t` e queda de pressão média `d_pressure` entre entrada e saída. | Todo passo |
| `flow_error.dat` | Duas colunas: tempo `t` e erro relativo percentual da vazão (`(Q_out - Q_in)/Q_in * 100`). | Todo passo |
| `area_error.dat` | Duas colunas: tempo `t` e erro relativo percentual da área da gota (`(A - A_ini)/A_ini * 100`). | Todo passo |
| `wall_contour.dat` | Coordenadas `(x, y)` dos nós do canal. | Passo 1 e a cada 10 passos |
| `drop_t=<t_ms>.dat` | Coordenadas `(x, y)` do contorno da gota no instante `nint(t * 1000)` em milésimos adimensionais (ex.: `drop_t=  125.dat` para `t = 0.125`). | Passo 1 e a cada 10 passos |
