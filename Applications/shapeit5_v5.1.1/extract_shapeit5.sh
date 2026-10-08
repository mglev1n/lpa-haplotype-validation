#!/bin/bash
# Copy the ShapeIt5 phase_common_static binary out of an existing container
# image into this folder, so the Docker build can install it from the
# repository. GitHub disabled odelaneau/shapeit5, so the binary can no longer be
# downloaded from its release page.
#
# Usage (from the repository root):
#   bash Applications/shapeit5_v5.1.1/extract_shapeit5.sh containers/lpa-validation-singularity_latest.sif

set -euo pipefail

SIF=${1:?Usage: $0 path/to/container.sif}
DEST_DIR="$(cd "$(dirname "$0")" && pwd)"
DEST="$DEST_DIR/phase_common_static"
SOURCE_PATH=/usr/local/bin/phase_common_static

[[ -f "$SIF" ]] || { echo "ERROR: container not found: $SIF" >&2; exit 1; }

if ! command -v apptainer &> /dev/null; then
    module load apptainer
fi

echo "Copying $SOURCE_PATH from $SIF"
apptainer exec --cleanenv "$SIF" cat "$SOURCE_PATH" > "$DEST.tmp"

# Reject anything that is not a Linux executable (for example, the HTML page
# that newer builds installed in place of the binary)
if [[ "$(head -c 4 "$DEST.tmp" | tail -c 3)" != "ELF" ]]; then
    echo "ERROR: $SOURCE_PATH in $SIF is not an ELF executable. First 200 bytes:" >&2
    head -c 200 "$DEST.tmp" >&2; echo >&2
    rm -f "$DEST.tmp"
    exit 1
fi

mv "$DEST.tmp" "$DEST"
chmod +x "$DEST"

SIZE_MB=$(( $(stat -c %s "$DEST") / 1024 / 1024 ))
if (( SIZE_MB >= 95 )); then
    echo "WARNING: binary is ${SIZE_MB} MB; GitHub rejects files over 100 MB" >&2
fi

echo
echo "Binary banner:"
"$DEST" --help 2>&1 | head -5 || true

# Record where the binary came from
{
    echo "phase_common_static provenance"
    echo "Extracted: $(date -u +'%Y-%m-%d %H:%M:%S UTC') by ${USER:-unknown} on $(hostname)"
    echo "Source container: $(realpath "$SIF")"
    echo "Source container modified: $(stat -c %y "$SIF")"
    echo "Source path in container: $SOURCE_PATH"
    echo "Container version file:"
    apptainer exec --cleanenv "$SIF" cat /opt/lpa-pipeline/VERSION 2>/dev/null | sed 's/^/  /' || echo "  (none)"
    echo "SHA-256: $(sha256sum "$DEST" | cut -d' ' -f1)"
    echo "Size: $(stat -c %s "$DEST") bytes"
} > "$DEST_DIR/PROVENANCE.txt"

echo
cat "$DEST_DIR/PROVENANCE.txt"
echo
echo "Next, from the repository root:"
echo "  git add Applications/shapeit5_v5.1.1/phase_common_static Applications/shapeit5_v5.1.1/PROVENANCE.txt"
echo "  git commit -m \"Add ShapeIt5 v5.1.1 phase_common_static binary\""
echo "  git push"
