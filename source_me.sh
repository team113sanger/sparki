#!/bin/bash

OPT_DIRECTORY="/lustre/scratch126/casm/team113da/users/jb62/opt"

export R_LIBS_USER="/lustre/scratch126/casm/team113da/users/jb62/rlibs"
export RENV_PATHS_ROOT="${OPT_DIRECTORY:?}/renv"
export RENV_PATHS_LIBRARY_ROOT="${OPT_DIRECTORY:?}/renv/library"
export RENV_PATHS_LIBRARY="${OPT_DIRECTORY:?}/renv/library/sparki-libs"
export RENV_PATHS_CACHE="${OPT_DIRECTORY:?}/renv-cache"

export PKGTYPE="binary"
export RENV_CONFIG_PAK_ENABLED="TRUE"

mkdir -p ${R_LIBS_USER}
mkdir -p ${RENV_PATHS_ROOT}
mkdir -p ${RENV_PATHS_LIBRARY}
mkdir -p ${RENV_PATHS_CACHE}

module load rocker/rver/4.4.0
module load git
module load pre-commit
