const core = require('@actions/core');
const github = require('@actions/github');
const exec = require('@actions/exec');

try {
    var filter = core.getInput('filter');
    if (!filter) {
        const ref = core.getInput('ref', {required: true});
        const prefix = core.getInput('prefix', {required: true});
        filter = `.ref == "${ref}" and (.key | startswith("${prefix}"))`;
    }
    const ids = await exec.getExecOutput("gh", [
        "cache", "list", "--sort", "created_at", "--json", "createdAt,id,key,ref", "--jq",
        `.[] | select(${filter}) | .id`
    ]).stdout.trim().split();
    const keep = core.getInput('keep');
    console.log(`Cache now contains ${ids.length} entries, deleting ${max(ids.length - keep, 0)} and keeping ${min(ids.length, keep)}`);
    for (const key of ids.slice(keep)) {
        console.log(`Deleting ${key}`);
        if (await exec.exec("gh", ["cache", "delete", key]) != 0) {
            core.error(`Could not delete cache entry with key ${key}`);
        }
    }
} catch (error) {
    core.setFailed(error.message);
}
