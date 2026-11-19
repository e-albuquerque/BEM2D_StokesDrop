# ==============================================================================
# Makefile - 2D Stokes Flow BEM Drop Simulation (Processor / Gota)
# Compilador: GNU Fortran (gfortran) + LAPACK/BLAS + OpenMP
# Detecta automaticamente Windows (Windows_NT) ou Linux/macOS
# ==============================================================================

FC      := gfortran
OPENMP  ?= 1

FFLAGS  := -O3 -march=native -fdefault-real-8 -fdefault-double-8 -ffree-line-length-512
LDFLAGS :=
LIBS    := -llapack -lblas

ifeq ($(OPENMP),1)
  FFLAGS  += -fopenmp
  LDFLAGS += -fopenmp
endif

SRC_DIR   := Processor
BUILD_DIR := build

ifeq ($(OS),Windows_NT)
  SHELL       := cmd.exe
  TARGET      := $(SRC_DIR)/Gota.exe
  LDFLAGS     += -static -Wl,--stack,536870912
  MKDIR_BUILD = if not exist "$(subst /,\,$(BUILD_DIR))" mkdir "$(subst /,\,$(BUILD_DIR))"
  RM_RF       = if exist "$(subst /,\,$(BUILD_DIR))" rmdir /s /q "$(subst /,\,$(BUILD_DIR))"
  RM_BIN      = if exist "$(subst /,\,$(TARGET))" del /q "$(subst /,\,$(TARGET))"
else
  TARGET      := $(SRC_DIR)/Gota
  MKDIR_BUILD = mkdir -p $(BUILD_DIR)
  RM_RF       = rm -rf $(BUILD_DIR)
  RM_BIN      = rm -f $(TARGET)
endif

MOD_OBJS := \
	$(BUILD_DIR)/variables.o \
	$(BUILD_DIR)/variables_build_HGB.o

SUB_OBJS := \
	$(BUILD_DIR)/alloc_memory.o \
	$(BUILD_DIR)/alloc_memory_build_HGB.o \
	$(BUILD_DIR)/apply_bc.o \
	$(BUILD_DIR)/apply_bc2.o \
	$(BUILD_DIR)/building_HGB.o \
	$(BUILD_DIR)/building_HGB_sing.o \
	$(BUILD_DIR)/calc_circle_curv.o \
	$(BUILD_DIR)/calc_curv_nodes.o \
	$(BUILD_DIR)/calc_drop_area.o \
	$(BUILD_DIR)/calc_radius.o \
	$(BUILD_DIR)/Conj_grad.o \
	$(BUILD_DIR)/lapack_solver.o \
	$(BUILD_DIR)/curv_sign.o \
	$(BUILD_DIR)/driver_pressure.o \
	$(BUILD_DIR)/evolution_euler.o \
	$(BUILD_DIR)/evolution_euler_implicit.o \
	$(BUILD_DIR)/first_drop_elem.o \
	$(BUILD_DIR)/flow_error.o \
	$(BUILD_DIR)/gauss_elimination.o \
	$(BUILD_DIR)/Gauss_Legendre.o \
	$(BUILD_DIR)/Gauss_Logaritimic.o \
	$(BUILD_DIR)/geometric_parameters.o \
	$(BUILD_DIR)/mat_prod_2x2_2x6.o \
	$(BUILD_DIR)/nodes_relocation.o \
	$(BUILD_DIR)/prod_Axb.o \
	$(BUILD_DIR)/read_input.o \
	$(BUILD_DIR)/regular_integral.o \
	$(BUILD_DIR)/reorganize.o \
	$(BUILD_DIR)/shape_functions.o \
	$(BUILD_DIR)/shape_functions_sing.o \
	$(BUILD_DIR)/sing_integral.o \
	$(BUILD_DIR)/sing_reg_integral.o \
	$(BUILD_DIR)/tecplot_contour.o \
	$(BUILD_DIR)/transform_eta_qsi.o \
	$(BUILD_DIR)/write_HGB.o \
	$(BUILD_DIR)/write_matA.o \
	$(BUILD_DIR)/Processing.o

OBJS := $(MOD_OBJS) $(SUB_OBJS)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90
	@$(MKDIR_BUILD)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

# Dependência explícita: todas as subrotinas aguardam a geração dos módulos (.mod)
$(SUB_OBJS): $(MOD_OBJS)

clean:
	$(RM_RF)
	$(RM_BIN)
