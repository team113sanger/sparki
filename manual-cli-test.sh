#!/bin/bash
set -eou pipefail

OUTDIR="manual-outdir"
rm -rf ${OUTDIR}
mkdir -p ${OUTDIR}

Rscript -e "devtools::load_all(); SPARKI::cli()"\
  --std-reports 'tests/testthat/testdata/std_reports/case_with_usual_reports' \
  --mpa-reports 'tests/testthat/testdata/mpa_reports/case_with_usual_reports' \
  --reference 'tests/testthat/testdata/inspect.txt' \
  --organism 'Homo sapiens' \
  --domain 'Viruses,Bacteria' \
  --outdir "${OUTDIR}" \
  --verbose
