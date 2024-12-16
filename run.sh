#!/usr/bin/env bash

snakemake --use-conda -p --rerun-incomplete --profile debug #--rerun-trigger mtime

