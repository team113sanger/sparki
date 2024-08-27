#!/bin/bash

_print_error() {
    local message="${1}"
    # The log message should include the name of the file that is currently
    # being executed
    echo -e "\033[0;31mERROR: ${message}\033[0m" >&2
}

_assert_file_exists() {
    local file_path="${1}"
    if [ ! -f "${file_path}" ]; then
        _print_error "File '${file_path}' does not exist"
        exit 1
    fi
}
