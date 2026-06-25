#!/bin/bash
# Double-click this in Finder to set up Locus. It just runs setup.sh in this folder.
exec "$(cd "$(dirname "$0")" && pwd)/setup.sh"
