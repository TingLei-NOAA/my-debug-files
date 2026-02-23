#!/usr/bin/env bash
set -euo pipefail

# ----------------------------
# HARD-WIRED SETTINGS
# ----------------------------
BASE_COMMIT="34354155"

# >>> EDIT THIS to your colleague-tree absolute path <<<
CHANGEDIR="/lfs/h2/emc/da/noscrub/Ting.Lei/dr-rdasapp/RDASApp/sorc/_workaround_/fv3-jedi-io"

# Destination root is where you run the script (should be .../fv3-jedi/src/fv3jedi)
DEST_ROOT="$(pwd)"

# ----------------------------
# Sanity checks
# ----------------------------
if [[ ! -d "$CHANGEDIR" ]]; then
  echo "ERROR: CHANGEDIR does not exist or is not a directory: $CHANGEDIR"
  exit 2
fi

# Need to be inside a git repo so git show BASE:PATH works
if ! git rev-parse --show-toplevel >/dev/null 2>&1; then
  echo "ERROR: Not inside a git repository (needed for git show BASE:PATH)."
  exit 2
fi

REPO_TOP="$(git rev-parse --show-toplevel)"

echo "============================================================"
echo "3-way merge colleague tree into DEST_ROOT"
echo "Repo top    : $REPO_TOP"
echo "Base commit : $BASE_COMMIT"
echo "Change dir  : $CHANGEDIR"
echo "Dest root   : $DEST_ROOT"
echo "============================================================"

if ! git cat-file -e "${BASE_COMMIT}^{commit}" >/dev/null 2>&1; then
  echo "ERROR: Base commit not found: $BASE_COMMIT"
  exit 2
fi

# DEST_ROOT relative to repo top (used to locate BASE files in git)
case "$DEST_ROOT" in
  "$REPO_TOP"/*) DEST_REL_TO_REPO="${DEST_ROOT#${REPO_TOP}/}" ;;
  *)
    echo "ERROR: DEST_ROOT is not inside the git repo top."
    echo "  DEST_ROOT=$DEST_ROOT"
    echo "  REPO_TOP=$REPO_TOP"
    exit 2
    ;;
esac

echo "Dest root relative to repo: $DEST_REL_TO_REPO"

# ----------------------------
# Main loop over files in CHANGEDIR
# ----------------------------
find "$CHANGEDIR" -type f -print0 | while IFS= read -r -d '' SRCFILE; do
  # Path relative to CHANGEDIR
  RELPATH="${SRCFILE#${CHANGEDIR}/}"

  # Destination file under DEST_ROOT
  DSTFILE="${DEST_ROOT}/${RELPATH}"

  # Repo-relative path to the file (needed for git show BASE_COMMIT:path)
  GITPATH="${DEST_REL_TO_REPO}/${RELPATH}"

  echo
  echo "------------------------------------------------------------"
  echo "Processing: $RELPATH"
  echo "  Source (THEIRS): $SRCFILE"
  echo "  Target (OURS)  : $DSTFILE"
  echo "  Base path      : $BASE_COMMIT:$GITPATH"

  # Ensure destination directory exists
  DSTDIRE="$(dirname "$DSTFILE")"
  if [[ ! -d "$DSTDIRE" ]]; then
    echo "  Creating directory: $DSTDIRE"
    mkdir -p "$DSTDIRE"
  fi

  # If target doesn't exist: copy directly
  if [[ ! -f "$DSTFILE" ]]; then
    echo "  Target missing -> copying new file"
    cp -p "$SRCFILE" "$DSTFILE"
    continue
  fi

  # Backup existing target
  BACKUP="${DSTFILE}-org"
  echo "  Backing up current file to: $BACKUP"
  cp -p "$DSTFILE" "$BACKUP"

  # If base commit doesn't contain that path, can't 3-way merge -> replace (backup kept)
  if ! git cat-file -e "${BASE_COMMIT}:${GITPATH}" >/dev/null 2>&1; then
    echo "  WARNING: Base commit lacks $GITPATH"
    echo "  -> Replacing with colleague version (backup kept)."
    cp -p "$SRCFILE" "$DSTFILE"
    continue
  fi

  TMPBASE="$(mktemp)"
  TMPTHEIRS="$(mktemp)"
  TMPMERGED="$(mktemp)"

  # BASE from git
  git show "${BASE_COMMIT}:${GITPATH}" > "$TMPBASE"
  # THEIRS from colleague tree
  cp -p "$SRCFILE" "$TMPTHEIRS"

  echo "  3-way merge (OURS=current, BASE=$BASE_COMMIT, THEIRS=colleague)..."

  if git merge-file -p "$DSTFILE" "$TMPBASE" "$TMPTHEIRS" > "$TMPMERGED"; then
    echo "  Merge completed."
  else
    echo "  WARNING: merge-file reported conflicts/issues; conflict markers may be present."
  fi

  cp -p "$TMPMERGED" "$DSTFILE"

  if grep -qE '^(<<<<<<<|=======|>>>>>>>)' "$DSTFILE"; then
    echo "  !!! CONFLICT MARKERS detected in $RELPATH"
    echo "  Resolve manually. Backup is: ${DSTFILE}-org"
  else
    echo "  No conflict markers detected."
  fi

  rm -f "$TMPBASE" "$TMPTHEIRS" "$TMPMERGED"
done

echo
echo "============================================================"
echo "Done. Next:"
echo "  cd \"$DEST_ROOT\""
echo "  git diff"
echo "  grep -R \"^<<<<<<<\" -n ."
echo "  git add -A && git commit -m \"Integrate colleague changes (3-way vs $BASE_COMMIT)\""
echo "============================================================"
