#############################################################
# Compile options start                                     #
#-----------------------------------------------------------#
#                                                           #

#----------------------------------------------#
# Library PATH and compiler settings           #
#----------------------------------------------#

CUDA_INSTALL_PATH ?= /usr/local/cuda
CUDA_SAMPLES_PATH ?= /usr/local/cuda/samples
FFTW_INSTALL_PATH ?= /usr/local/fftw-3.3.8
CPPCOMPILER       ?= g++
MPICOMPILER       ?= mpicxx
OPTIMIZATION      ?= -O3
OMPFLAG           ?= -fopenmp
 # (If you use g++ compiler, please set the value as "-fopenmp".)

#----------------------------------------------#
# Environment settings                         #
#----------------------------------------------#

USE_GPU    := 1
 # (If you do not use GPU, please set the value as 0 or comment out of the line.)

USE_MPI    := 0
 # (If you do not use MPI, please set the value as 0 or comment out of the line.) 

#                                                           #
#-----------------------------------------------------------#
# Compile options end                                       #
#############################################################

#-----------------------------------------------------------------------------------------#
CUFILES    := fft_process.cu main.cu
CCFILES    := control.cpp cpu_time.cpp docking.cpp parameter.cpp pdb_entry.cpp \
              protein.cpp receptor.cpp ligand.cpp

ifndef USE_GPU
       USE_GPU := 0
endif
ifndef USE_MPI
       USE_MPI := 0
endif

ifeq ($(USE_GPU),1)
     CXXFLAGS   += -DCUFFT
     NVCCFLAGS  += -DCUFFT
     ifeq ($(USE_MPI),1)
          CCFILES    += mpidp.cpp
          EXECUTABLE := megadock-gpu-dp
          COMPILER    = $(MPICOMPILER)
          CXXFLAGS   += -DMPI_DP
          NVCCFLAGS  += -DMPI_DP
          OBJDIRPR   := gm
     else     
          EXECUTABLE := megadock-gpu
          COMPILER    = $(CPPCOMPILER)
          OBJDIRPR   := gs
     endif
else
     CUFILES :=
     CCFILES += fft_process.cpp main.cpp
     ifeq ($(USE_MPI),1)
          CCFILES    += mpidp.cpp
          EXECUTABLE := megadock-dp
          COMPILER    = $(MPICOMPILER)
          CXXFLAGS   += -DMPI_DP
          NVCCFLAGS  += -DMPI_DP
          OBJDIRPR   := cm
     else
          EXECUTABLE := megadock
          COMPILER    = $(CPPCOMPILER)
          OBJDIRPR   := cs
     endif
endif

.SUFFIXES : .cuda.c .cu_dbg.o .c_dbg.o .cpp_dbg.o .cu_rel.o .c_rel.o .cpp_rel.o .cubin .ptx

# FFTW setup
FFTW_CFLAGS  ?= -I$(FFTW_INSTALL_PATH)/include
FFTW_LDFLAGS ?= -lm -L$(FFTW_INSTALL_PATH)/lib -lfftw3f

# Add new SM Versions here as devices with new Compute Capability are released
SM_VERSIONS := sm_60

# Basic directory setup for SDK
SRCDIR     ?= src
ROOTDIR    ?= $(CUDA_SAMPLES_PATH)
BINDIR     ?= .
ROOTOBJDIR ?= obj_$(OBJDIRPR)
COMMONDIR  := $(ROOTDIR)/common

# Compilers
NVCC       := $(CUDA_INSTALL_PATH)/bin/nvcc -Xcompiler -fopenmp -arch=$(SM_VERSIONS) -use_fast_math
CXX        := $(COMPILER) $(OMPFLAG)
LINK       := $(COMPILER) -fPIC $(OMPFLAG)

# Includes
INCLUDES  += -I$(SRCDIR)
ifeq ($(USE_GPU),1)
	INCLUDES += -I$(CUDA_INSTALL_PATH)/include -I$(COMMONDIR)/inc
endif

# Warning flags
CXXWARN_FLAGS := #	-W -Wall \
#	-Wimplicit \
#	-Wformat \
#	-Wparentheses \
#	-Wmultichar \
#	-Wtrigraphs \
#	-Wpointer-arith \
#	-Wreturn-type \
#	-Wno-unused-function

# Compiler-specific flags
NVCCFLAGS += 
CXXFLAGS  += -static $(CXXWARN_FLAGS)

# Common flags
COMMONFLAGS += $(OPTIMIZATION) $(INCLUDES) $(FFTW_CFLAGS)
LIBSUFFIX   := _x86_64
NVCCFLAGS   += --compiler-options -fno-strict-aliasing
#CXXFLAGS    += -fno-strict-aliasing -funroll-loops

# FFTW Libs
LIB       := $(FFTW_LDFLAGS)

# CUDA Libs
ifeq ($(USE_GPU),1)
	LIB   += -L$(CUDA_INSTALL_PATH)/lib64 -lcudart -lcufft
endif

TARGETDIR := $(BINDIR)
TARGET    := $(TARGETDIR)/$(EXECUTABLE)
LINKLINE   = $(LINK) -o $(TARGET) $(OBJS) $(LIB)

################################################################################
# Check for input flags and set compiler flags appropriately
################################################################################
ifdef maxregisters
	NVCCFLAGS += -maxrregcount $(maxregisters)
endif

# Add cudacc flags
NVCCFLAGS += $(CUDACCFLAGS)

# Add common flags
NVCCFLAGS += $(COMMONFLAGS)
CXXFLAGS  += $(COMMONFLAGS)

################################################################################
# Set up object files
################################################################################
OBJDIR := $(ROOTOBJDIR)
OBJS +=  $(patsubst %.cpp,$(OBJDIR)/%.cpp.o,$(notdir $(addprefix $(SRCDIR)/, $(CCFILES))))
OBJS +=  $(patsubst %.c,$(OBJDIR)/%.c.o,$(notdir $(addprefix $(SRCDIR)/, $(CCFILES))))
OBJS +=  $(patsubst %.cu,$(OBJDIR)/%.cu.o,$(notdir $(addprefix $(SRCDIR)/, $(CUFILES))))
OBJS :=  $(filter %.o,$(OBJS))

################################################################################
# Rules
################################################################################

$(OBJDIR)/%.c.o : $(SRCDIR)/%.c $(C_DEPS)
	$(VERBOSE)$(CXX) $(CXXFLAGS) -o $@ -c $<

$(OBJDIR)/%.cpp.o : $(SRCDIR)/%.cpp $(C_DEPS)
	$(VERBOSE)$(CXX) $(CXXFLAGS) -o $@ -c $<

$(OBJDIR)/%.cu.o : $(SRCDIR)/%.cu $(CU_DEPS)
	$(VERBOSE)$(NVCC) $(NVCCFLAGS) $(SMVERSIONFLAGS) -o $@ -c $<

define SMVERSION_template
OBJS += $(patsubst %.cu,$(OBJDIR)/%.cu_$(1).o,$(notdir $(CUFILES_$(1))))
$(OBJDIR)/%.cu_$(1).o : $(SRCDIR)%.cu $(CU_DEPS)
	$(VERBOSE)$(NVCC) -o $$@ -c $$< $(NVCCFLAGS) -arch $(1)
endef

# This line invokes the above template for each arch version stored in
# SM_VERSIONS.  The call funtion invokes the template, and the eval
# function interprets it as make commands.
$(foreach smver,$(SM_VERSIONS),$(eval $(call SMVERSION_template,$(smver))))

$(TARGET): messeages makedirectories $(OBJS) decoygen calcrg Makefile
	$(VERBOSE)$(LINKLINE)

messeages :
	@echo --------------------------
	@echo make $(EXECUTABLE) start
	@echo --------------------------

makedirectories :
	$(VERBOSE)mkdir -p $(OBJDIR)

decoygen : $(SRCDIR)/decoygen.cpp
	$(VERBOSE)$(CPPCOMPILER) $(SRCDIR)/decoygen.cpp -lm -o decoygen

calcrg : $(SRCDIR)/calcrg.cpp
	$(VERBOSE)$(CPPCOMPILER) $(SRCDIR)/calcrg.cpp -o calcrg

.PHONY : clean allclean
clean :
	$(VERBOSE)rm -rf $(ROOTOBJDIR)
	$(VERBOSE)rm -f $(TARGET)

allclean :
	$(VERBOSE)rm -rf obj_cs obj_gs obj_cm obj_gm
	$(VERBOSE)rm -f megadock megadock-gpu megadock-dp megadock-gpu-dp decoygen calcrg
