#!/usr/bin/env python3
import os
import re
import argparse
import sys
import subprocess

try:
    import yaml
except ImportError:
    print("Error: PyYAML is required to run this script. Please install it using: pip install PyYAML")
    sys.exit(1)

def load_yaml(path):
    if not os.path.exists(path):
        return {}
    with open(path) as fh:
        return yaml.safe_load(fh) or {}

def replace_vars(d, context):
    if isinstance(d, dict):
        return {k: replace_vars(v, context) for k, v in d.items()}
    elif isinstance(d, list):
        return [replace_vars(x, context) for x in d]
    elif isinstance(d, str):
        val = d
        # Loop to resolve recursive templates
        for _ in range(10):
            matches = re.findall(r'\{\{([A-Za-z0-9_]+)\}\}', val)
            if not matches:
                break
            resolved_any = False
            for m in matches:
                if m in context:
                    val = val.replace("{{" + m + "}}", str(context[m]))
                    resolved_any = True
            if not resolved_any:
                break
        return val
    return d

def find_containers(d, path_prefix=""):
    containers = []
    if isinstance(d, dict):
        for k, v in d.items():
            current_path = f"{path_prefix}.{k}" if path_prefix else k
            if k in ["container", "default_container"] and isinstance(v, str):
                containers.append((current_path, v))
            else:
                containers.extend(find_containers(v, current_path))
    elif isinstance(d, list):
        for i, v in enumerate(d):
            containers.extend(find_containers(v, f"{path_prefix}[{i}]"))
    return containers

def main():
    parser = argparse.ArgumentParser(description="Resolve and download Singularity/Apptainer containers for offline use.")
    parser.add_argument("--config", default="config/config.yaml", help="Path to main config.yaml")
    parser.add_argument("--marvin-config", default="config/site_configs/site_config_marvin.yaml", help="Path to Marvin site config (to resolve Docker URIs)")
    parser.add_argument("--miarka-config", default="config/site_configs/site_config_miarka.yaml", help="Path to Miarka site config (to resolve local SIF names)")
    parser.add_argument("--out-dir", help="Override output directory for SIF files (default: APPTAINER_CACHE from Miarka config)")
    parser.add_argument("--dry-run", action="store_true", help="Only show the download commands without executing them")
    parser.add_argument("--pull-tool", choices=["apptainer", "singularity"], default="apptainer", help="Tool to use for pulling containers")
    args = parser.parse_args()

    # Load configurations
    main_config = load_yaml(args.config)
    marvin_site = load_yaml(args.marvin_config)
    miarka_site = load_yaml(args.miarka_config)

    if not main_config:
        print(f"Error: Could not load main config from {args.config}")
        sys.exit(1)

    # Merge contexts for interpolation
    marvin_context = {**main_config, **marvin_site}
    miarka_context = {**main_config, **miarka_site}

    # Interpolate
    resolved_marvin = replace_vars(main_config, marvin_context)
    resolved_miarka = replace_vars(main_config, miarka_context)

    # Find containers
    marvin_containers = find_containers(resolved_marvin)
    miarka_containers = find_containers(resolved_miarka)

    container_mappings = {}
    for (path_marvin, uri), (path_miarka, sif_path) in zip(marvin_containers, miarka_containers):
        if path_marvin != path_miarka:
            continue
        # Skip if it doesn't look like a template or docker URI
        if not uri.startswith("docker://"):
            continue
        
        # Determine output filename
        sif_filename = os.path.basename(sif_path)
        if not sif_filename.endswith(".sif"):
            # Try to construct standard name if not resolved to SIF extension
            sif_filename = sif_filename.replace(":", "_").replace("/", "_") + ".sif"
            
        container_mappings[uri] = sif_filename

    # Determine apptainer cache output dir
    out_dir = args.out_dir or miarka_site.get("APPTAINER_CACHE")
    if not out_dir:
        # Fallback if no APPTAINER_CACHE defined
        out_dir = "apptainer_cache"

    print(f"Resolving containers from: {args.config}")
    print(f"Output directory: {out_dir}\n")

    if not os.path.exists(out_dir) and not args.dry_run:
        os.makedirs(out_dir, exist_ok=True)

    for uri, filename in sorted(container_mappings.items()):
        dest_path = os.path.join(out_dir, filename)
        cmd = [args.pull_tool, "pull", dest_path, uri]
        print(f"Pulling: {uri} -> {dest_path}")
        print(f"Command: {' '.join(cmd)}")
        
        if not args.dry_run:
            if os.path.exists(dest_path):
                print(f"File already exists: {dest_path}, skipping.\n")
                continue
            try:
                subprocess.run(cmd, check=True)
                print("Successfully pulled!\n")
            except subprocess.CalledProcessError as e:
                print(f"Error pulling container {uri}: {e}\n")
        else:
            print()

if __name__ == "__main__":
    main()
