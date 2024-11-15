# syntax=docker/dockerfile:1
# MUTLI-STAGE DOCKERFILE for building a Python image
# STAGE 1: base_stage is the base image for both the development and dev-staging/staging/production images
# STAGE 2: development_only_stage is the image used for development (optimised for VSCode)

########################
# STAGE 1: base_stage #
########################

# IMPORTANT
# If you change the base image, you will need to update the
# PRE_FETCH_BASE_IMAGE variable in the .gitlab-ci.yml file also.
FROM rocker/r-ver:4.4.0 AS base_stage

USER root
# Arguments for UID and GID matching
ARG USER_ID=1000
ARG GROUP_ID=1000

# Set the top level environment variables
# PKGTYPE - prefer renv/instal.packages/usethis::use_package to prefer binary packages
# RENV_CONFIG_PAK_ENABLED - enable pak as the middleware for package management
ENV \
    DATA_DIRECTORY="/data" \
    OPT_DIRECTORY="/opt" \
    USER_NAME="admin" \
    USER_DIRECTORY="/home/admin" \
    LC_ALL="en_US.UTF-8" \
    LANG="en_US.UTF-8" \
    PKGTYPE="binary" \
    RENV_CONFIG_PAK_ENABLED="TRUE" \
    PRE_COMMIT_VERSION="4.0.1"


# Set next environment variables that interpolate the top level environment
# variables
#
# To build isolated R environments for packaging projects we want R deps to be
# installed outside of the project directory (this also helps with bind
# mounting). This principally effects the RENV_PATHS_LIBRARY_ROOT. For more
# information see https://rstudio.github.io/renv/articles/packages.html#library-paths
ENV \
    USER_BASHRC="${USER_DIRECTORY:?}/.bashrc" \
    USER_BIN_DIRECTORY="${USER_DIRECTORY:?}/.local/bin" \
    SSH_DIR="${USER_DIRECTORY:?}/.ssh" \
    PROJECT_DIRECTORY="${OPT_DIRECTORY:?}/repo" \
    RENV_DIRECTORY="${OPT_DIRECTORY:?}/renv" \
    RENV_PATHS_ROOT="${OPT_DIRECTORY:?}/renv" \
    RENV_PATHS_LIBRARY_ROOT="${OPT_DIRECTORY:?}/renv/library" \
    RENV_PATHS_LIBRARY="${OPT_DIRECTORY:?}/renv/library/sparki-libs" \
    RENV_PATHS_CACHE="${OPT_DIRECTORY:?}/renv-cache" \
    PIPX_HOME="${OPT_DIRECTORY}/pipx" \
    PIPX_BIN_DIR="${OPT_DIRECTORY}/pipx/bin" \
    LOGGING_DIRECTORY="${DATA_DIRECTORY:?}/logs"

# Set the environment that interpolates the next level of environment variables
#
# Set the PATH to include the pipx bin directory
ENV \
    PATH=${PIPX_BIN_DIR:?}:${PATH}

# Run the commands to:
# - create directories defined in the environment variables
# - create the docker group if it does not exist (helpful for bind mounting)
# - create the non-root user
# - put the non-root user in the same groups as the docker user
# - give non-root user ownership of the directories
RUN \
    locale-gen "${LANG:?}" \
    && update-locale LANG="${LANG:?}" \
    && groupadd -g ${GROUP_ID} "${USER_NAME}" \
    && useradd -u ${USER_ID} -g ${GROUP_ID} "${USER_NAME}" --shell /bin/bash --create-home --home-dir "${USER_DIRECTORY}" \
    && if ! getent group docker > /dev/null; then groupadd docker; fi \
    && usermod -a -G docker,staff admin \
    && mkdir -p "${PROJECT_DIRECTORY:?}" "${DATA_DIRECTORY:?}" "${OPT_DIRECTORY:?}" "${RENV_DIRECTORY:?}" "${RENV_PATHS_LIBRARY:?}" "${RENV_PATHS_CACHE:?}" "${PIPX_BIN_DIR:?}" "${PIPX_BIN_DIR:?}" \
    && chown -R "${USER_NAME:?}:${USER_NAME:?}" "${PROJECT_DIRECTORY:?}" "${DATA_DIRECTORY:?}" "${USER_DIRECTORY:?}" "${OPT_DIRECTORY:?}" \
    && chmod -R 755 "${PROJECT_DIRECTORY:?}" "${DATA_DIRECTORY:?}" "${USER_DIRECTORY:?}" "${OPT_DIRECTORY:?}"

################################################
# Install the system dependencies              #
################################################
#
# For compiling and installing the htslib, samtools, libdeflate, ruby we need the following packages:
# - curl, wget, tree, nano for general development
# - vim (specifically vi) is needed by R usethis
# - build-essential for compiling
# - git for version control
# - pipx to install python built tools, pre-commit (apt-get version of pre-commit is too old)
# themselves -- this is done in a later stage.
RUN \
    apt-get update -y \
    && DEBIAN_FRONTEND=noninteractive apt-get install -yq --no-install-recommends \
    build-essential \
    curl \
    wget \
    tree \
    vim \
    nano \
    git \
    pipx \
    && pipx install "pre-commit==${PRE_COMMIT_VERSION:?}" \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

####################################################
# Install R packages for bootstrapping development #
#                                                  #
# Use Renv+pak via a wrapper script                #
####################################################

# We use the renv package to bootstrap the R package installation and isolate them from the system
# R installation.
RUN R -e 'install.packages("renv")'

# Copy the renv.lock file (if it exists) and the DESCRIPTION file to the
# container (as well as the renv/ directory)
WORKDIR $PROJECT_DIRECTORY
COPY --chown="${USER_NAME}:${USER_NAME}" ["DESCRIPTION", "NAMESPACE", "./"]

# Install the R packages from the lock file or the DESCRIPTION file if the lock
# file does not exist.
RUN if [ -f renv.lock ]; then \
    Rscript -e 'renv::init(bare = TRUE); renv::restore()' ; \
    else \
    Rscript -e 'renv::init(bare = TRUE); renv::install()' ; \
    fi
# Copy the rest of the project files to the container
COPY --chown="${USER_NAME}:${USER_NAME}" . .

# Reapply permissions after all installation and copying is done so the user can
# manipulate the files if necessary
RUN \
    chown -R "${USER_NAME:?}:${USER_NAME:?}" "${PROJECT_DIRECTORY:?}" "${DATA_DIRECTORY:?}" "${USER_DIRECTORY:?}" "${OPT_DIRECTORY:?}" \
    && chmod -R 755 "${PROJECT_DIRECTORY:?}" "${DATA_DIRECTORY:?}" "${USER_DIRECTORY:?}" "${OPT_DIRECTORY:?}"

USER "${USER_NAME:?}"


###################################
# STAGE 2: development_only_stage #
# - This stage is optional        #
# - It is optimised for VSCode    #
###################################

# To develop from the container we need add some extra directories to work nicely with VSCode
FROM base_stage AS development_only_stage
USER root

# Install hubflow to allow for git flow style development & conditional install
# sudo, giving the user passwordless sudo privileges
WORKDIR "${USER_DIRECTORY}"
ARG HAS_SUDO="${HAS_SUDO:-0}"
RUN git config --global --add safe.directory "${USER_DIRECTORY}/gitflow" \
    && git config --global --add safe.directory "${USER_DIRECTORY}/gitflow/shFlags" \
    && git clone https://github.com/datasift/gitflow \
    && chown -R "${USER_NAME}:${USER_NAME}" gitflow \
    && cd gitflow \
    && ./install.sh \
    && if [ "${HAS_SUDO}" = "1" ]; then \
    apt-get update -y \
    && apt-get install -y sudo \
    && rm -rf /var/lib/apt/lists/* \
    && echo "${USER_NAME:?} ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers; \
    fi

# Return back to the project directory for the entrypoint
WORKDIR "${PROJECT_DIRECTORY}"
USER "${USER_NAME}"

# Prepare the directory for the VSCode server extensions and bash history
RUN mkdir -p "${USER_DIRECTORY}/.vscode-server/extensions" \
    "${USER_DIRECTORY}/.vscode-server-insiders/extensions" \
    && chown -R "${USER_NAME}:${USER_NAME}" \
    "${USER_DIRECTORY}/.vscode-server" \
    "${USER_DIRECTORY}/.vscode-server-insiders" && \
    SNIPPET="export PROMPT_COMMAND='history -a' && export HISTFILE=${USER_DIRECTORY}/.commandhistory/.bash_history" \
    && mkdir "${USER_DIRECTORY}/.commandhistory/" \
    && touch "${USER_DIRECTORY}/.commandhistory/.bash_history" \
    && chown -R "${USER_NAME}:${USER_NAME}" "${USER_DIRECTORY}/.commandhistory/" \
    && echo "$SNIPPET" >> "/home/$USER_NAME/.bashrc"

# Install pre-commit hooks
RUN pre-commit install --install-hooks
