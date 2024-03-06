const core = require('@actions/core');
const github = require('@actions/github');
const exec = require('@actions/exec');

// https://github.com/actions/cache/blob/main/tips-and-workarounds.md
// https://docs.github.com/en/actions/using-workflows/caching-dependencies-to-speed-up-workflows#force-deleting-cache-entries

try {
    var filter = core.getInput('filter');
    if (!filter) {
        const ref = core.getInput('ref');
        const prefix = core.getInput('prefix');
        filter = `.ref == "${ref}" and (.key | startswith("${prefix}"))`;
    }
    const ids = exec.getExecOutput("gh", [
        "cache", "list", "--sort", "created_at", "--json", "createdAt,id,key,ref", "--jq",
        `.[] | select(${filter}) | .id`
    ]).stdout.trim().split();
    console.log(`Cache contains ${ids.length} old entries`);
} catch (error) {
    core.setFailed(error.message);
}
