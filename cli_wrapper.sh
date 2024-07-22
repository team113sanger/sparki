#!/bin/bash

set -eu

# CONSTANTS

ENV_VAR__FUNCNAME="CLI_FUNCTION_NAME"
ENV_VAR__PROG_NAME="PROG_NAME"
# We need to find the cli.R file which is in the same directory as this script
CODE_DIRECTORY="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")"
DERMATLAS_R_BIN=/software/team113/dermatlas/R/R-4.2.2/bin
CLI_FILE="$(readlink -f "${CODE_DIRECTORY}/src/cli.R")"

# FUNCTIONS

# A function to print in red and to stderr an error message
_print_error() {
    local message="${1}"
    # The log message should include the name of the file that is currently
    # being executed
    local current_file="${BASH_SOURCE[1]##*/}"
    echo -e "\033[0;31mERROR ($current_file): ${message}\033[0m" >&2
}
# A function to check if a command exists
_assert_command_exists() {
    local command_name="${1}"
    if ! command -v "${command_name}" &> /dev/null; then
        _print_error "Command '${command_name}' does not exist"
        exit 1
    fi
}

# A function to check if a directory exists and is readable
_assert_dir_exists() {
    local dir_path="${1}"
    if [ ! -d "${dir_path}" ]; then
        _print_error "Directory '${dir_path}' does not exist"
        exit 1
    fi

    if [ ! -r "${dir_path}" ]; then
        _print_error "Directory '${dir_path}' is not readable"
        exit 1
    fi
}

# A function to check if a file exists
_assert_file_exists() {
    local file_path="${1}"
    if [ ! -f "${file_path}" ]; then
        _print_error "File '${file_path}' does not exist"
        exit 1
    fi
}

# A function to check if an environment variable is set
_assert_env_var_set() {
    local env_var_name="${1}"
    if [ -z "${!env_var_name:-}" ]; then
        _print_error "Environment variable '${env_var_name}' is not set"
        exit 1
    fi
}

# RUNTIME

_assert_dir_exists "${CODE_DIRECTORY}"
_assert_dir_exists "${DERMATLAS_R_BIN}"
_assert_file_exists "${CLI_FILE}"
_assert_env_var_set "${ENV_VAR__FUNCNAME}"
_assert_env_var_set "${ENV_VAR__PROG_NAME}"
export PATH="${DERMATLAS_R_BIN}:${PATH}" 
_assert_command_exists "Rscript"

Rscript "${CLI_FILE}" "${!ENV_VAR__FUNCNAME}" "$@"
