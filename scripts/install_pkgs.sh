#!/bin/bash
# Installs the Python requirements at the start of Claude Code cloud sessions.
# On your own computer this exits right away; there, run
#   python -m pip install -r requirements.txt
# once yourself.

if [ "$CLAUDE_CODE_REMOTE" != "true" ]; then
  exit 0
fi

cd "$CLAUDE_PROJECT_DIR" || exit 0
log=/tmp/install_pkgs.log

# pip skips requirements that are already satisfied. The second command covers
# images whose system Python refuses plain pip installs (PEP 668).
if pip install -q -r requirements.txt >"$log" 2>&1 \
   || pip install -q --break-system-packages -r requirements.txt >>"$log" 2>&1; then
  echo "SessionStart: Python requirements are installed."
else
  echo "SessionStart: pip install failed (see $log). Install requirements.txt by hand."
fi
exit 0
