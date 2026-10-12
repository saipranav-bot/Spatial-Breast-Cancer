#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

OUT="data/raw/outs"
H5="${OUT}/filtered_feature_bc_matrix.h5"
SP="${OUT}/spatial"

mkdir -p "${OUT}" "${SP}"

if [[ -e "${H5}" ]] && [[ -d "${SP}" ]] && [[ -n "$(ls -A "${SP}" 2>/dev/null)" ]]; then
  echo "Data already present in ${OUT}; skipping download."
  exit 0
fi

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

echo "Downloading 10x Visium V1_Breast_Cancer_Block_A_Section_1..."
curl -L -o "${TMP}/V1_Breast_Cancer_Block_A_Section_1_raw_feature_bc_matrix.h5" \
  "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Breast_Cancer_Block_A_Section_1/V1_Breast_Cancer_Block_A_Section_1_raw_feature_bc_matrix.h5"
curl -L -o "${TMP}/V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5" \
  "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Breast_Cancer_Block_A_Section_1/V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5"
curl -L -o "${TMP}/V1_Breast_Cancer_Block_A_Section_1_spatial.tar.gz" \
  "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Breast_Cancer_Block_A_Section_1/V1_Breast_Cancer_Block_A_Section_1_spatial.tar.gz"

tar -xzf "${TMP}/V1_Breast_Cancer_Block_A_Section_1_spatial.tar.gz" -C "${TMP}"

cp "${TMP}/V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5" "${H5}"
cp "${TMP}/spatial/"* "${SP}/" 2>/dev/null || true

echo "Downloaded to ${OUT}"
