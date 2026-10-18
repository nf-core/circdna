#!/usr/bin/env bash
# prepare_eccsplorer_database.sh
#
# Merge Dfam and RepBase annotation databases into a single FASTA file
# for ECCsplorer BLAST annotation.
#
# Usage:
#   bash scripts/prepare_eccsplorer_database.sh --species fungi
#   bash scripts/prepare_eccsplorer_database.sh --species human --out /path/to/output.fa
#
# Species mapping (RepBase *.ref files):
#   fungi|fungus|yeast -> fngrep.ref
#   plant|arabidopsis|ath -> athrep.ref
#   human|hum|homo -> humrep.ref
#   <other> -> ${species}rep.ref
set -euo pipefail

# -----------------------------------------------------------------------------
# Default configuration
# -----------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

DFAM_DIR="/Users/siyangming/Library/CloudStorage/OneDrive-个人/Project_backup/PublicDB/Dfam"
DFAM_GZ="${DFAM_DIR}/Dfam-RepeatMasker.lib.gz"
DFAM_LIB="${DFAM_DIR}/Dfam-RepeatMasker.lib"

REPBASE_DIR="/Users/siyangming/Library/CloudStorage/OneDrive-个人/Project_backup/PublicDB/RepBase/RepBase31.07.fasta"

OUT_DIR="${PIPELINE_DIR}/testdatasets/eccsplorer_db"
OUT_FILE=""

SPECIES=""

# -----------------------------------------------------------------------------
# Parse command-line arguments
# -----------------------------------------------------------------------------
usage() {
    cat <<EOF
Usage: $(basename "$0") --species <name> [--out <path>]

Options:
    --species <name>   Species name for RepBase selection.
                    Mappings:
                        fungi|fungus|yeast -> fngrep.ref
                        plant|arabidopsis|ath -> athrep.ref
                        human|hum|homo -> humrep.ref
                        <other> -> \${species}rep.ref
    --out <path>       Output file or directory path.
                    If a directory, writes eccsplorer_db.fa inside it.
                    Default: ${OUT_DIR}/eccsplorer_db.fa
    -h, --help          Show this help message.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --species)
            SPECIES="$2"
            shift 2
            ;;
        --out)
            OUT_FILE="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "[ERROR] Unknown argument: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ -z "${SPECIES}" ]]; then
    echo "[ERROR] --species is required." >&2
    usage >&2
    exit 1
fi

# -----------------------------------------------------------------------------
# Resolve species -> RepBase filename
# -----------------------------------------------------------------------------
case "${SPECIES}" in
    fungi|fungus|yeast)
        REP_FILE="fngrep.ref"
        ;;
    plant|arabidopsis|ath)
        REP_FILE="athrep.ref"
        ;;
    human|hum|homo)
        REP_FILE="humrep.ref"
        ;;
    *)
        REP_FILE="${SPECIES}rep.ref"
        ;;
esac

REP_PATH="${REPBASE_DIR}/${REP_FILE}"

# -----------------------------------------------------------------------------
# Resolve output path
# -----------------------------------------------------------------------------
if [[ -z "${OUT_FILE}" ]]; then
    OUT_FILE="${OUT_DIR}/eccsplorer_db.fa"
elif [[ -d "${OUT_FILE}" ]]; then
    OUT_FILE="${OUT_FILE}/eccsplorer_db.fa"
fi

OUT_FINAL_DIR="$(dirname "${OUT_FILE}")"
mkdir -p "${OUT_FINAL_DIR}"

# -----------------------------------------------------------------------------
# Validate inputs
# -----------------------------------------------------------------------------
echo "[INFO] Species:      ${SPECIES}"
echo "[INFO] RepBase file: ${REP_PATH}"
echo "[INFO] Output:       ${OUT_FILE}"

if [[ ! -f "${DFAM_GZ}" && ! -f "${DFAM_LIB}" ]]; then
    echo "[ERROR] Dfam database not found." >&2
    echo "        Looked for: ${DFAM_LIB}" >&2
    echo "        Looked for: ${DFAM_GZ}" >&2
    exit 1
fi

if [[ ! -f "${REP_PATH}" ]]; then
    echo "[ERROR] RepBase file not found: ${REP_PATH}" >&2
    echo "        Available *.ref files in ${REPBASE_DIR}:" >&2
    ls -1 "${REPBASE_DIR}"/*.ref 2>/dev/null | sed 's/^/          /' >&2 || true
    exit 1
fi

# -----------------------------------------------------------------------------
# Decompress Dfam if needed
# -----------------------------------------------------------------------------
# Prefer a pre-decompressed copy in the original Dfam directory (read-only use).
# If it is not there (or the directory is read-only), decompress into the output
# directory so the cache is reusable across runs.
DFAM_LIB_LOCAL="${OUT_FINAL_DIR}/Dfam-RepeatMasker.lib"

if [[ -f "${DFAM_LIB}" ]]; then
    echo "[INFO] Dfam already decompressed (source dir): ${DFAM_LIB}"
    DFAM_USE="${DFAM_LIB}"
elif [[ -f "${DFAM_LIB_LOCAL}" ]]; then
    echo "[INFO] Dfam already decompressed (output dir): ${DFAM_LIB_LOCAL}"
    DFAM_USE="${DFAM_LIB_LOCAL}"
else
    echo "[INFO] Decompressing Dfam database..."
    echo "        ${DFAM_GZ}"
    echo "        -> ${DFAM_LIB_LOCAL}"
    gunzip -c "${DFAM_GZ}" > "${DFAM_LIB_LOCAL}"
    echo "[INFO] Decompression complete."
    DFAM_USE="${DFAM_LIB_LOCAL}"
fi

# -----------------------------------------------------------------------------
# Concatenate databases into a single FASTA file
# -----------------------------------------------------------------------------
echo "[INFO] Concatenating Dfam + RepBase (${REP_FILE})..."
cat "${DFAM_USE}" "${REP_PATH}" > "${OUT_FILE}"

# -----------------------------------------------------------------------------
# Verify output
# -----------------------------------------------------------------------------
if [[ ! -s "${OUT_FILE}" ]]; then
    echo "[ERROR] Output file is empty: ${OUT_FILE}" >&2
    exit 1
fi

SEQ_COUNT=$(grep -c '^>' "${OUT_FILE}" || true)
echo "[INFO] Output written: ${OUT_FILE}"
echo "[INFO] File size:      $(du -h "${OUT_FILE}" | cut -f1)"
echo "[INFO] Sequence count: ${SEQ_COUNT}"
echo "[INFO] First 3 FASTA headers:"
awk '/^>/{print; c++} c==3{exit}' "${OUT_FILE}" | sed 's/^/        /'
echo "[INFO] Done."
