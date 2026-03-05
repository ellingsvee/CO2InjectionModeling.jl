#!/bin/bash
set -e
cd "$(dirname "$0")/.."
docker run --rm \
  -v "$(pwd)/joss:/data" \
  -e JOURNAL=joss \
  openjournals/inara \
  -o pdf paper.md
echo "Output: joss/paper.pdf"
