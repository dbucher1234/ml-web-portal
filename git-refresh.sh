#!/usr/bin/env bash
#
# git-refresh.sh  —  one-command add → commit → pull-rebase → push
# Usage: ./git-refresh.sh "my commit message"
# If no argument is given, the commit message defaults to "Quick sync".

set -e  # abort on first error

msg=${1:-"Quick sync"}                    # use arg1 or default message
branch=$(git symbolic-ref --short HEAD)   # current branch name (master, main, dev…)

echo "🔄  Syncing branch: $branch"

git add -A
# commit quietly if nothing to commit
if ! git diff --cached --quiet; then
    git commit -m "$msg"
else
    echo "🛈  Nothing new to commit."
fi

# ensure we have an upstream; set it the first time if missing
if ! git rev-parse --abbrev-ref --symbolic-full-name "@{u}" &>/dev/null; then
    echo "⤴️  Setting upstream origin/$branch"
    git push -u origin "$branch"
else
    git pull  origin "$branch" --rebase
    git push  origin "$branch"
fi

echo "✅  Done."

