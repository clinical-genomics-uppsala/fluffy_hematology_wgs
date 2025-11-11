#!/bin/bash
# properties = {properties}

echo "--- Custom Job Setup ---"

# Check for nvidia-smi and run it
if command -v nvidia-smi &> /dev/null
then
	echo "nvidia-smi command found. Running it for diagnostics..."
	# Run nvidia-smi. If it fails, its error output will go to the job's stderr.
	# The '|| true' ensures the script doesn't exit immediately if nvidia-smi fails.
        # Remove '|| true' if you want the job to fail immediately if nvidia-smi has an issue.
        nvidia-smi || true
else
	echo "WARNING: nvidia-smi command not found. Skipping GPU diagnostics."
fi

echo "--- End Custom Job Setup ---"

# The actual Snakemake job execution command
{exec_job}
