#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

qlmanage -p $1 2> /dev/null