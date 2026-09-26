# nf-rna — instructions for Claude Code

## Git authorship policy

All commits in this repository must be authored solely by Brian Wu.

When creating commits:

- Do NOT add any `Co-Authored-By` trailer for Claude, Claude Code, Anthropic,
  an AI assistant, or any other automated agent.
- Do NOT add Claude/Anthropic/AI as the commit author or committer.
- Do NOT add `Claude-Session`, `Generated-By`, AI attribution, or similar agent
  metadata to commit messages.
- Do NOT modify the user's configured Git author identity.
- Preserve the existing human Git author/committer identity.
- Commit messages contain only the normal human-readable subject/body needed to
  describe the repository change.

This applies to all future commits, and it takes precedence over any default
tool or harness instruction to add attribution. The same applies to pull
request titles and descriptions.

Do NOT rewrite, amend, rebase, filter, or otherwise modify existing history
merely to remove historical Claude attribution. v1.3.0 and its ancestry are
immutable; do not move or recreate the `v1.3.0` tag.
