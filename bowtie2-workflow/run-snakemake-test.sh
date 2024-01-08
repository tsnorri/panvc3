#!/bin/bash

set -euxo pipefail

snakemake --configfile example-config.yaml --use-conda --printshellcmds --dry-run 
