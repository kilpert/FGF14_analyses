#!/usr/bin/env bash

snakemake --cores 32 --use-conda -p --rerun-incomplete --profile humgen #--rerun-trigger mtime

