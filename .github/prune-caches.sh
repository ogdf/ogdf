#!/bin/bash

caches=$(gh cache list --key "$1" --ref "$GITHUB_REF" --repo "$GITHUB_REPOSITORY" --sort created_at --order desc --json id --jq ".[].id")
count=$(echo "$caches" | grep -v ^$ | wc -l)
if [ "$count" -eq 0 ]; then
    echo "Found no caches for $GITHUB_REPOSITORY:$GITHUB_REF:$1 to prune."
elif [ "$count" -eq 1 ]; then
    echo "Found 1 cache entry ($caches) for $GITHUB_REPOSITORY:$GITHUB_REF:$1, keeping it."
else
    echo "Found $count cache entries for $GITHUB_REPOSITORY:$GITHUB_REF:$1:"
    echo $caches
    echo "Only retaining the latest (first) one."
    echo "$caches" | tail -n +2 | xargs -L 1 -P 4 --no-run-if-empty --verbose -exec gh cache delete
fi
