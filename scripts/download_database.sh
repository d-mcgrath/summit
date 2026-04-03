#!/usr/bin/env bash

set -euo pipefail

DOI_URL="https://doi.org/10.5281/zenodo.19390972"
RECORD_ID="19390973"
ARCHIVE_NAME="merged_regulatory_features_bgz.tar.gz"
DEST_ROOT="."
FORCE="false"

usage() {
  cat <<'EOF'
Download and unpack the Summit merged_regulatory_features bundle from Zenodo.

Usage:
  bash scripts/download_database.sh [options]

Options:
  --doi-url URL            Zenodo DOI URL.
  --record-id ID           Zenodo record ID.
  --archive-name NAME      Archive filename stored on Zenodo.
  --dest-root PATH         Destination root directory. Default: current directory.
  --force true|false       Re-download and overwrite existing extracted files. Default: false.
  --help                   Show this help text.

Defaults are configured for:
  DOI:        https://doi.org/10.5281/zenodo.19390972
  Record ID:  19390973
  Archive:    merged_regulatory_features_bgz.tar.gz

Typical usage from the repository root:
  bash scripts/download_database.sh --dest-root .

After extraction, Summit expects:
  data/merged_regulatory_features_bgz/
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --doi-url)
      DOI_URL="$2"
      shift 2
      ;;
    --record-id)
      RECORD_ID="$2"
      shift 2
      ;;
    --archive-name)
      ARCHIVE_NAME="$2"
      shift 2
      ;;
    --dest-root)
      DEST_ROOT="$2"
      shift 2
      ;;
    --force)
      FORCE="$2"
      shift 2
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

if ! command -v curl >/dev/null 2>&1; then
  echo "curl is required but was not found on PATH." >&2
  exit 1
fi

if ! command -v tar >/dev/null 2>&1; then
  echo "tar is required but was not found on PATH." >&2
  exit 1
fi

discover_archive_from_api() {
  local api_url payload
  api_url="https://zenodo.org/api/records/${RECORD_ID}"

  if command -v python3 >/dev/null 2>&1; then
    payload="$(curl -L --fail --silent --show-error "$api_url" | python3 - "$ARCHIVE_NAME" <<'PY'
import json, sys
preferred = sys.argv[1]
try:
    data = json.load(sys.stdin)
except Exception:
    sys.exit(1)
files = data.get("files", []) or []
if not files:
    sys.exit(1)

def pick(files):
    for item in files:
        key = item.get("key", "")
        if key == preferred:
            return item
    for item in files:
        key = item.get("key", "")
        if key.endswith(".tar.gz") and "merged_regulatory_features" in key:
            return item
    for item in files:
        key = item.get("key", "")
        if key.endswith(".tar.gz"):
            return item
    for item in files:
        key = item.get("key", "")
        if key.endswith(".zip") or key.endswith(".tgz"):
            return item
    return None

item = pick(files)
if item is None:
    sys.exit(1)
links = item.get("links", {}) or {}
url = links.get("self") or links.get("content") or ""
if not url:
    sys.exit(1)
print(item.get("key", ""))
print(url)
PY
)" || return 1
  elif command -v python >/dev/null 2>&1; then
    payload="$(curl -L --fail --silent --show-error "$api_url" | python - "$ARCHIVE_NAME" <<'PY'
import json, sys
preferred = sys.argv[1]
try:
    data = json.load(sys.stdin)
except Exception:
    sys.exit(1)
files = data.get("files", []) or []
if not files:
    sys.exit(1)

def pick(files):
    for item in files:
        key = item.get("key", "")
        if key == preferred:
            return item
    for item in files:
        key = item.get("key", "")
        if key.endswith(".tar.gz") and "merged_regulatory_features" in key:
            return item
    for item in files:
        key = item.get("key", "")
        if key.endswith(".tar.gz"):
            return item
    for item in files:
        key = item.get("key", "")
        if key.endswith(".zip") or key.endswith(".tgz"):
            return item
    return None

item = pick(files)
if item is None:
    sys.exit(1)
links = item.get("links", {}) or {}
url = links.get("self") or links.get("content") or ""
if not url:
    sys.exit(1)
print(item.get("key", ""))
print(url)
PY
)" || return 1
  else
    return 1
  fi

  printf '%s' "$payload"
}

discover_archive_name() {
  local record_url html candidates preferred
  record_url="https://zenodo.org/records/${RECORD_ID}"

  html="$(curl -L --fail --silent --show-error "$record_url")" || return 1

  candidates="$(printf '%s' "$html" | grep -oE "/records/${RECORD_ID}/files/[^\"?]+" | sed "s#^/records/${RECORD_ID}/files/##" | sort -u || true)"

  if [[ -z "$candidates" ]]; then
    return 1
  fi

  preferred="$(printf '%s\n' "$candidates" | grep -E 'merged_regulatory_features.*\.tar\.gz$' | head -n 1 || true)"
  if [[ -n "$preferred" ]]; then
    printf '%s' "$preferred"
    return 0
  fi

  preferred="$(printf '%s\n' "$candidates" | grep -E '\.tar\.gz$' | head -n 1 || true)"
  if [[ -n "$preferred" ]]; then
    printf '%s' "$preferred"
    return 0
  fi

  preferred="$(printf '%s\n' "$candidates" | grep -E '\.(zip|tgz)$' | head -n 1 || true)"
  if [[ -n "$preferred" ]]; then
    printf '%s' "$preferred"
    return 0
  fi

  return 1
}

DEST_ROOT="$(cd "$DEST_ROOT" && pwd)"
TARGET_DIR="$DEST_ROOT/data/merged_regulatory_features_bgz"

if [[ -d "$TARGET_DIR" && "$FORCE" != "true" ]]; then
  echo "Merged regulatory features bundle already exists at:"
  echo "  $TARGET_DIR"
  echo "Use --force true to re-download and overwrite."
  exit 0
fi

DOWNLOAD_URL="https://zenodo.org/records/${RECORD_ID}/files/${ARCHIVE_NAME}?download=1"
TMP_DIR="$(mktemp -d)"
ARCHIVE_PATH="$TMP_DIR/$ARCHIVE_NAME"

cleanup() {
  rm -rf "$TMP_DIR"
}
trap cleanup EXIT

echo "Downloading merged regulatory features bundle"
echo "  DOI:      $DOI_URL"
echo "  Record:   $RECORD_ID"
echo "  Archive:  $ARCHIVE_NAME"
echo "  URL:      $DOWNLOAD_URL"
echo "  Dest:     $DEST_ROOT"

if ! curl -L --fail --output "$ARCHIVE_PATH" "$DOWNLOAD_URL"; then
  echo "Initial download attempt failed for:"
  echo "  $DOWNLOAD_URL"
  echo "Trying to discover the archive link from the Zenodo API..."

  API_DISCOVERY="$(discover_archive_from_api || true)"

  if [[ -n "${API_DISCOVERY:-}" ]]; then
    ARCHIVE_NAME="$(printf '%s\n' "$API_DISCOVERY" | sed -n '1p')"
    DOWNLOAD_URL="$(printf '%s\n' "$API_DISCOVERY" | sed -n '2p')"
    ARCHIVE_PATH="$TMP_DIR/$ARCHIVE_NAME"

    echo "Discovered archive via API:"
    echo "  $ARCHIVE_NAME"
    echo "Retrying download..."
    echo "  $DOWNLOAD_URL"

    curl -L --fail --output "$ARCHIVE_PATH" "$DOWNLOAD_URL"
  else
    echo "Zenodo API discovery failed. Trying to discover the archive filename from the record page..."

    DISCOVERED_ARCHIVE_NAME="$(discover_archive_name || true)"

    if [[ -z "${DISCOVERED_ARCHIVE_NAME:-}" ]]; then
      echo "Could not auto-detect an archive filename from the Zenodo record." >&2
      echo "Please check the record page and rerun with:" >&2
      echo "  --archive-name <exact_uploaded_filename>" >&2
      exit 1
    fi

    ARCHIVE_NAME="$DISCOVERED_ARCHIVE_NAME"
    ARCHIVE_PATH="$TMP_DIR/$ARCHIVE_NAME"
    DOWNLOAD_URL="https://zenodo.org/records/${RECORD_ID}/files/${ARCHIVE_NAME}?download=1"

    echo "Discovered archive name:"
    echo "  $ARCHIVE_NAME"
    echo "Retrying download..."
    echo "  $DOWNLOAD_URL"

    curl -L --fail --output "$ARCHIVE_PATH" "$DOWNLOAD_URL"
  fi
fi

mkdir -p "$DEST_ROOT"

if [[ "$FORCE" == "true" && -d "$TARGET_DIR" ]]; then
  echo "Removing existing extracted directory:"
  echo "  $TARGET_DIR"
  rm -rf "$TARGET_DIR"
fi

echo "Extracting archive..."
tar -xzf "$ARCHIVE_PATH" -C "$DEST_ROOT"

if [[ -d "$TARGET_DIR" ]]; then
  :
elif [[ -d "$DEST_ROOT/merged_regulatory_features_bgz" ]]; then
  mkdir -p "$DEST_ROOT/data"
  if [[ -d "$TARGET_DIR" ]]; then
    rm -rf "$TARGET_DIR"
  fi
  mv "$DEST_ROOT/merged_regulatory_features_bgz" "$TARGET_DIR"
else
  echo "Extraction completed, but the expected directory was not found." >&2
  echo "Looked for:" >&2
  echo "  $TARGET_DIR" >&2
  echo "  $DEST_ROOT/merged_regulatory_features_bgz" >&2
  exit 1
fi

FILE_COUNT="$(find "$TARGET_DIR" -maxdepth 1 -type f | wc -l | tr -d ' ')"
BGZ_COUNT="$(find "$TARGET_DIR" -maxdepth 1 -type f -name '*.tsv.bgz' | wc -l | tr -d ' ')"
TBI_COUNT="$(find "$TARGET_DIR" -maxdepth 1 -type f -name '*.tsv.bgz.tbi' | wc -l | tr -d ' ')"

echo "Download and extraction complete."
echo "  Installed directory: $TARGET_DIR"
echo "  Total files:         $FILE_COUNT"
echo "  .bgz files:          $BGZ_COUNT"
echo "  .tbi files:          $TBI_COUNT"

if [[ "$BGZ_COUNT" -eq 0 || "$TBI_COUNT" -eq 0 ]]; then
  echo "Warning: expected bgzip/tabix files were not detected after extraction." >&2
  exit 1
fi

echo
echo "You can now run Summit with:"
echo "  Rscript scripts/run_summit.R --config=config/example_merged_regulatory_features.yml"
