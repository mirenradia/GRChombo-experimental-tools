# -*- Mode: Makefile -*-

# the location of the Chombo "lib" directory
ifndef CHOMBO_HOME
    $(error Please define CHOMBO_HOME - see installation instructions.)
endif

# GRChombo Source directory
GRCHOMBO_SOURCE = ../../../GRChombo/Source

# trace the chain of included makefiles
makefiles += releasedExamples_AMRPoisson_execVariableCoefficient

## Define the variables needed by Make.example

# the base name(s) of the application(s) in this directory
ebase = Main_PoissonSolver

# names of Chombo libraries needed by this program, in order of search.
LibNames = AMRElliptic AMRTools BoxTools

# the locations of the source code directories
base_dir = .

# input file for 'run' target
INPUT = params.txt

# application-specific targets
src_dirs := Source \
						$(GRCHOMBO_SOURCE)/utils \
						$(GRCHOMBO_SOURCE)/InitialConditions/BosonStars/SingleStar

.DEFAULT_GOAL := all-git

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

src_dirs_full := $(PWD) $(src_dirs)
git_status_check := ""

git_status_check += $(foreach src_dir,$(src_dirs_full),$(shell git diff --stat $(src_dir)/*.[c,h]pp))
ifneq ($(strip $(git_status_check)),"")
	git_status := "-dirty"
endif
git_status ?= ""

git_commit_hash := $(shell git rev-parse --short HEAD)

all-git: all
	$(ECHO)mv $(_app_configs) $(ebase)$(config).$(git_commit_hash)$(git_status).ex
