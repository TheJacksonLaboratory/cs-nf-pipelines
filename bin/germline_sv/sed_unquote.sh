#!/usr/bin/env bash

# This is a very simple script to just remove double-quotes from files.
# It's broken out as a script to avoid issues with unescaped quotes
# interferring with the Nextflow script block.

sed -i 's/\"//g;s/\;//g' ${1}