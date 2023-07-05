#!/bin/bash

set -euxo pipefail

ulimit -n 4096
exec "$@"
