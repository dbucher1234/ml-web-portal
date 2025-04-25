#!/bin/bash
git add .
git commit -m "Quick sync" 2>/dev/null
git pull origin master --rebase   # <- change main → master
git push origin master            # <- change main → master

