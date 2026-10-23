# Copy this file, adapt CONDA_BASE, and export its absolute path before running any run_*.sh:
#   cp envs/env_setup.example.sh ~/dartseq_env.sh;  export ENV_SETUP=~/dartseq_env.sh
# Every job script then runs  . "$ENV_SETUP" <tools|r|bullseye>  instead of `module load ...`.
# Without ENV_SETUP the scripts fall back to the cluster's `module load` lines (written for the cluster they were developed on).
CONDA_BASE=${CONDA_BASE:-$HOME/miniforge3}          # where conda/mamba is installed
. "$CONDA_BASE/etc/profile.d/conda.sh" || return 1
case "$1" in
  tools)    conda activate dartseq-tools ;;
  r)        conda activate dartseq-r ;;
  bullseye) conda activate dartseq-bullseye ;;
  *) echo "env_setup: unknown environment '$1' (tools | r | bullseye)" >&2; return 1 ;;
esac
