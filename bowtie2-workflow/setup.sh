#!/bin/bash

set -euxo pipefail

mkdir -p input
gunzip -c -k ../test-workflow/merged.vcf.gz > input/merged.test.vcf
gunzip -c -k ../test-workflow/genome/genome.fa.gz > input/genome.fa
