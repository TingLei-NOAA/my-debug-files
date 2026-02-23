#!/usr/bin/env bash
set -euo pipefail

# Revert all permission/mode changes back to what HEAD records,
# while keeping ALL content changes (staged or unstaged).
# Prints warnings if chmod is not permitted.

git rev-parse --is-inside-work-tree >/dev/null 2>&1 || {
  echo "Error: not inside a git repository." >&2
  exit 1
}

head_mode() {
  local p="$1"
  git ls-tree -r HEAD -- "$p" | awk 'NR==1{print $1}'
}

fix_one() {
  local p="$1"
  local hm perm

  hm="$(head_mode "$p")"
  if [[ -z "$hm" ]]; then
    # Not in HEAD (new file). Skip.
    return 0
  fi

  perm="${hm:3:3}"   # 100644 -> 644, 100755 -> 755

  # 1) Fix index executable bit (Git tracks only +x/-x)
  if [[ "$hm" == "100755" ]]; then
    git update-index --chmod=+x -- "$p" >/dev/null 2>&1 || true
  else
    git update-index --chmod=-x -- "$p" >/dev/null 2>&1 || true
  fi

  # 2) Fix working tree file permissions (if file exists)
  if [[ -e "$p" ]]; then
    if ! chmod "$perm" -- "$p"; then
      echo "WARNING: chmod $perm failed for $p (permissions/mount/ACL may prevent chmod)" >&2
    fi
  fi
}

# Collect mode-changed paths from unstaged + staged summaries.
# Lines look like:
#  mode change 100644 => 100755 path
collect_mode_paths() {
  git "$@" --summary |
    awk '
      $1=="mode" && $2=="change" {
        # path is last field; for safety, print from field 6 onward if any spaces (rare)
        # Standard format: mode change <old> => <new> <path>
        for (i=6; i<=NF; i++) {
          printf "%s%s", $i, (i<NF ? OFS : "\n")
        }
      }
    '
}

tmp="$(mktemp)"
trap 'rm -f "$tmp"' EXIT

{
  collect_mode_paths diff
  collect_mode_paths diff --cached
} | sort -u > "$tmp"

# Apply fixes
while IFS= read -r p; do
  [[ -n "$p" ]] || continue
  fix_one "$p"
done < "$tmp"

echo "Remaining mode diffs (unstaged):"
git diff --summary | grep "mode change" || echo "  (none)"

echo "Remaining mode diffs (staged):"
git diff --cached --summary | grep "mode change" || echo "  (none)"

