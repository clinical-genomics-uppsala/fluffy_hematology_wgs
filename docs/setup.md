# Setup

This repository now follows the shared WP2 layout:

- `workflow/` for the Snakefile, rules, scripts, and schemas
- `config/` for pipeline config, references, resources, and outputs
- `profiles/local/` for workstation or MacBook runs
- `profiles/slurm/` for generic cluster submission

## Local run on MacBook M5

1. Install Snakemake in a dedicated environment, for example with `mamba` or `pipx`.
2. Install `hydra-genetics` in the same environment.
3. Install a container runtime. On macOS that usually means Docker Desktop plus a local Apptainer or Singularity-compatible setup if you want to execute containerized rules.
4. Clone the repository and point the workflow at local module checkouts with `module_root` if you do not want network access during module resolution.
5. Use `--profile profiles/local` for dry-runs or local execution.

## Offline packaging

For an offline cluster install, keep these pieces together:

- local copies of Hydra module repositories
- reference bundles under `config/references/`
- cluster-specific resources under `config/resources/`
- container images cached locally and referenced from the profile

The repository still supports the existing site-specific config files during the transition.
