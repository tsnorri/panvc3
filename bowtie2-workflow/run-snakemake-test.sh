#!/bin/bash

set -euxo pipefail

snakemake --configfile example-config.yaml --use-conda --conda-prefix ./conda-env --printshellcmds -c 16 --dry-run 
