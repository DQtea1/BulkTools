#!/usr/bin/env bash
# macOS launcher for the BulkTools Shiny app (equivalent of windows_launcher.bat).
# Double-click in Finder, or run "./macos_launcher.command" from a terminal.
set -euo pipefail

# --- Config -----------------------------------------------------------------
BROWSE_DIR="$HOME"            # Default root for selecting input files
OUT_DIR="$HOME/shiny_out"     # Where results are written (edit if you want)
IMAGE="qtea1/bulktools:latest"
PORT=5288

mkdir -p "$OUT_DIR"

# --- Allocate ~3/4 of the logical CPU cores to the container ---
CORES_TOTAL="$(sysctl -n hw.ncpu 2>/dev/null || echo 1)"
CORES=$(( CORES_TOTAL * 3 / 4 ))
[ "$CORES" -lt 1 ] && CORES=1
echo "Using ${CORES} of ${CORES_TOTAL} CPU core(s)."

# --- Make sure Docker Desktop is running ------------------------------------
if ! docker info >/dev/null 2>&1; then
  echo "Docker is not running. Trying to start Docker Desktop..."
  open -a Docker || { echo "Could not start Docker Desktop. Please open it manually, then re-run."; exit 1; }
  echo -n "Waiting for Docker to be ready"
  for _ in $(seq 1 60); do
    if docker info >/dev/null 2>&1; then echo " - OK"; break; fi
    echo -n "."
    sleep 1
  done
  if ! docker info >/dev/null 2>&1; then
    echo
    echo "Docker did not start in time. Open Docker Desktop and re-run this launcher."
    exit 1
  fi
fi

# --- Open the browser shortly after the app starts --------------------------
( sleep 4; open "http://localhost:${PORT}/" ) &

# External drives on macOS are mounted under /Volumes. Binding it makes every
# drive browsable in the app as /mounts/<DriveName> (see build_shiny_roots()).
VOLUMES_MOUNT=()
if [ -d /Volumes ]; then
  VOLUMES_MOUNT=(--mount type=bind,source=/Volumes,target=/mounts)
fi

# --- Run the app ------------------------------------------------------------
# --platform linux/amd64 forces Rosetta emulation on Apple Silicon (M1/M2/M3);
# it is a no-op on Intel Macs. The image is built for amd64 only.
docker run --rm -p ${PORT}:${PORT} \
  --platform linux/amd64 \
  --cpus=${CORES} \
  -e SHINY_PORT=${PORT} \
  -e SHINY_ROOT_PATH=/browse \
  -e SHINY_ROOT_NAME=home \
  -e SHINY_N_CORES=${CORES} \
  --mount type=bind,source="${BROWSE_DIR}",target=/browse \
  --mount type=bind,source="${OUT_DIR}",target=/out \
  "${VOLUMES_MOUNT[@]}" \
  "${IMAGE}"