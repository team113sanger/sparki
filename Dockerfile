# syntax=docker/dockerfile:1

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
    PKGTYPE="binary"

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
    LOGGING_DIRECTORY="${DATA_DIRECTORY:?}/logs"

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
    && mkdir -p "${PROJECT_DIRECTORY:?}" "${DATA_DIRECTORY:?}" "${OPT_DIRECTORY:?}" \
    && chown -R "${USER_NAME:?}:${USER_NAME:?}" "${PROJECT_DIRECTORY:?}" "${DATA_DIRECTORY:?}" "${USER_DIRECTORY:?}" "${OPT_DIRECTORY:?}" \
    && chmod -R 755 "${PROJECT_DIRECTORY:?}" "${DATA_DIRECTORY:?}" "${USER_DIRECTORY:?}" "${OPT_DIRECTORY:?}"

################################################
# Install the system dependencies              #
################################################
#
# For compiling and installing the htslib, samtools, libdeflate, ruby we need the following packages:
# - vim (specifically vi) is needed by R usethis
# - build-essential for compiling
# themselves -- this is done in a later stage.
RUN \
    apt-get update -y \
    && DEBIAN_FRONTEND=noninteractive apt-get install -yq --no-install-recommends \
    build-essential \
    vim \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

################################################
# Install R packages                           #
################################################

WORKDIR ${PROJECT_DIRECTORY}

# Copy the rest of the project files to the container
COPY --chown="${USER_NAME}:${USER_NAME}" . .

# We use the renv package to bootstrap the R package installation and isolate them from the system
# R installation.
RUN \
    rm -f .Rprofile \
    && rm -f renv.lock \
    && rm -rf renv/* \
    && R -e 'install.packages("devtools")' \
    && R -e 'devtools::install(dependencies = TRUE)'

# Reapply permissions after all installation and copying is done so the user can
# manipulate the files if necessary
RUN \
    chown -R "${USER_NAME:?}:${USER_NAME:?}" "${PROJECT_DIRECTORY:?}" "${DATA_DIRECTORY:?}" "${USER_DIRECTORY:?}" "${OPT_DIRECTORY:?}" \
    && chmod -R 755 "${PROJECT_DIRECTORY:?}" "${DATA_DIRECTORY:?}" "${USER_DIRECTORY:?}" "${OPT_DIRECTORY:?}"

USER "${USER_NAME:?}"
