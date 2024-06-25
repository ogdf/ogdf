#!/bin/bash

# use CLI args or default from env / working dir
KEY="$1"
REF="${2:-${WORKFLOW_REF:-${GITHUB_REF:-$(git symbolic-ref HEAD)}}}"
REPO="${3:-${GITHUB_REPOSITORY:-$(gh repo view --json nameWithOwner --jq .nameWithOwner)}}"

if [[ -z "$KEY" ]]; then
    echo "Need to specify a cache key as first arg!"
    exit 1
fi

caches=$(gh cache list --key "$KEY" --ref "$REF" --repo "$REPO" --sort created_at --order desc --json id --jq ".[].id")
count=$(echo "$caches" | grep -v ^$ | wc -l)
if [ "$count" -eq 0 ]; then
    echo "Found no caches for $REPO:$REF:$KEY to prune."
elif [ "$count" -eq 1 ]; then
    echo "Found 1 cache entry ($caches) for $REPO:$REF:$KEY, keeping it."
else
    echo "Found $count cache entries for $REPO:$REF:$KEY:"
    echo $caches
    echo "Only retaining the latest (first) one."
    echo "$caches" | tail -n +2 | xargs -L 1 -P 4 --no-run-if-empty --verbose -exec gh cache delete
fi
