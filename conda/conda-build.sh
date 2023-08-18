#!/bin/bash

set -euxo pipefail

conda build -c conda-forge .
