# AGENT.md

You are developing a production-quality bioinformatics software package.

General principles

- Python >=3.12
- Use type hints everywhere.
- Prefer dataclasses where appropriate.
- Every public function must have docstrings.
- Follow SOLID principles.
- Modular design.
- Avoid duplicate code.
- Use pathlib instead of os.path.
- Use logging instead of print().
- Write pytest tests for all important functions.
- Never write large monolithic scripts.

Libraries

Use

- pysam
- pyranges
- pandas
- numpy
- matplotlib
- gffutils
- tqdm
- loguru

Avoid shell commands whenever possible.

Performance

- BAM files may contain millions of reads.
- Do not load an entire BAM into memory.
- Use indexed interval queries.
- Minimize repeated overlap searches.

Coding style

- Black formatting
- Ruff linting
- Explicit imports
- Small functions
- Keep files under ~400 lines whenever possible.

Testing

Every new feature should include tests.

Do not consider work complete unless tests pass.

Version Control

This project uses Git.

Before modifying code:

1. Check git status.
2. Confirm current branch.

Rules:

- Do not commit directly to main.
- Feature-branch workflow: start each unit of work from `main` with
  `git switch -c <topic>` (e.g. `feat/step-07`, `fix/strand-overlap`),
  commit in small focused units, run tests before each commit, then merge
  back with `git switch main && git merge --no-ff <topic>`. Keep the topic
  branch after merging (do not delete it). There is no remote; merges are
  local only.
- Keep main clean and always-passing (tests green after every merge).
- Commit frequently, in small focused units (each commit passes tests).
- Run tests before committing.
- Use Conventional Commits for all commit messages:
  feat:, fix:, docs:, style:, refactor:, test:, chore:, perf:, build:, ci:.

License

- This project is released under the MIT License; a LICENSE file is required.

Python Environment Requirements

This project uses a dedicated Python virtual environment.

Before running Python commands:

- Use the project virtual environment.
- Do not install packages into the system Python environment.
- Do not modify global Python packages.

Expected environment:

.venv/
