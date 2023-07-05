#!/bin/bash

set -euxo pipefail

snakemake --printshellcmds --use-conda --conda-prefix ./conda-env --cores 16
