# syntax=docker/dockerfile:1

################################################
# Set up the system                            #
################################################

# IMPORTANT
# If you change the base image, you will need to update the
# PRE_FETCH_BASE_IMAGE variable in the .gitlab-ci.yml file also.
FROM rocker/r-ver:4.4.0 AS base_stage

USER root

#**** Set build & runtime variables ****#

ENV \
    OPT_DIRECTORY="/opt" \
    USER_NAME="admin" \
    USER_DIRECTORY="/home/admin" \
    LC_ALL="en_US.UTF-8" \
    LANG="en_US.UTF-8"

#**** Set build-time variables ****#

ARG USER_ID=1000
ARG GROUP_ID=1000
ARG PROJECT_DIRECTORY="${OPT_DIRECTORY:?}/repo"

# Run the commands to:
# - Create directories defined in the variables defined above;
# - Create the docker group if it does not exist (helpful for bind mounting);
# - Create the non-root user;
# - Put the non-root user in the same groups as the docker user;
# - Give non-root user ownership of the directories.
RUN \
    locale-gen "${LANG:?}" \
    && update-locale LANG="${LANG:?}" \
    && groupadd -g ${GROUP_ID} "${USER_NAME}" \
    && useradd -u ${USER_ID} -g ${GROUP_ID} "${USER_NAME}" --shell /bin/bash --create-home --home-dir "${USER_DIRECTORY}" \
    && if ! getent group docker > /dev/null; then groupadd docker; fi \
    && usermod -a -G docker,staff admin \
    && mkdir -p "${PROJECT_DIRECTORY:?}" "${OPT_DIRECTORY:?}" \
    && chown -R "${USER_NAME:?}:${USER_NAME:?}" "${PROJECT_DIRECTORY:?}" "${USER_DIRECTORY:?}" "${OPT_DIRECTORY:?}" \
    && chmod -R 755 "${PROJECT_DIRECTORY:?}" "${USER_DIRECTORY:?}" "${OPT_DIRECTORY:?}"

################################################
# Install the system dependencies              #
################################################

# Install a few programs in the system.
RUN \
    apt-get update -y \
    && DEBIAN_FRONTEND=noninteractive apt-get install -yq --no-install-recommends \
    build-essential \
    vim \
    nano \
    tree \
    curl \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

################################################
# Install SPARKI and its dependencies          #
################################################

WORKDIR ${PROJECT_DIRECTORY}

# Copy the files/directories necessary for SPARKI to be installed.
COPY --chown="${USER_NAME}:${USER_NAME}" "R" "R"
COPY --chown="${USER_NAME}:${USER_NAME}" "man" "man"
COPY --chown="${USER_NAME}:${USER_NAME}" "vignettes" "vignettes"
COPY --chown="${USER_NAME}:${USER_NAME}" ["NAMESPACE", "DESCRIPTION", "./"]

# Install devtools and then SPARKI from the local repository.
RUN R -e "install.packages('devtools'); devtools::install(dependencies = TRUE)"

# Reapply permissions after all installation and copying is done so the user can manipulate the files if necessary.
RUN \
    chown -R "${USER_NAME:?}:${USER_NAME:?}" "${PROJECT_DIRECTORY:?}" "${USER_DIRECTORY:?}" "${OPT_DIRECTORY:?}" \
    && chmod -R 755 "${PROJECT_DIRECTORY:?}" "${USER_DIRECTORY:?}" "${OPT_DIRECTORY:?}"

USER "${USER_NAME:?}"
