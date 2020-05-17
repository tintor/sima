#!/bin/bash
set -e

core/core_test -d=yes
geom/geom_test -d=yes
santorini/santorini_test -d=yes
./test -d=yes
