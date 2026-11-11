#!/usr/bin/env bash
# Waits on arm 1c's PID (not a name match -- a name match would catch this
# script's own pgrep and never terminate), then starts the KF9 arm.
set -uo pipefail
WATCH=3147
echo "waiting on 1c pid $WATCH, started $(date -u '+%F %T UTC')"
while [ "$WATCH" -ne 0 ] && kill -0 "$WATCH" 2>/dev/null; do sleep 60; done
echo "1c finished at $(date -u '+%F %T UTC'); starting KF9"
exec bash "$HOME/pmhc/run_kf9.sh"
