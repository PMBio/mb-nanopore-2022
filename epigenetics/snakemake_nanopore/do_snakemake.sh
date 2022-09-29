#!/bin/bash

if [[ -z "${SNAKEMAKE_PATH}" ]]; then
  snakemake="snakemake"
else
  snakemake="${SNAKEMAKE_PATH}"
fi

if [[ -z "${MEDUSA_CLUSTER_SCRIPT}" ]]; then
  echo "Please set environment variable MEDUSA_CLUSTER_SCRIPT." >&2;
  exit 1
else
  cluster_script="${MEDUSA_CLUSTER_SCRIPT}"
fi

if [[ -z "${SNAKEMAKE_CLUSTER_STATUS_SCRIPT}" ]]; then
  echo "SNAKEMAKE_CLUSTER_STATUS_SCRIPT not set. Jobs killed due to memory limits or timeouts will stall snakemake. Look at the subdirectory 'status' for example cluster status scripts."
  cluster_status_script=""
else
  cluster_status_script="--cluster-status ${SNAKEMAKE_CLUSTER_STATUS_SCRIPT}"
fi

if [[ -z "${HDF5_PLUGIN_PATH}" ]]; then
  echo "Please set environment variable HDF5_PLUGIN_PATH to a path that contains the vbz plugin. See https://github.com/nanoporetech/vbz_compression." >&2;
  exit 1
fi

${snakemake} --cluster "${cluster_script}" ${cluster_status_script} --jobs 128 --latency-wait 120 $@

