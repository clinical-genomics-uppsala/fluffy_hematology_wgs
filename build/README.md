
```bash
TAG_OR_BRANCH="v2.5.0" #adjust version 
PIPELINE_NAME="fluffy_hematology_wgs" 
PIPELINE_GITHUB_REPO="https://github.com/clinical-genomics-uppsala/fluffy_hematology_wgs"
PYTHON_VERSION="3.9"
bash build/build_conda.sh config/references/<REFERENCES.YAML> #adjust to correct file(s)
```
