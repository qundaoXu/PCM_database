#!/bin/zsh
set -e
pcm_project_dir="${0:A:h}"
exec /opt/anaconda3/envs/pcm_database/bin/python "$pcm_project_dir/scripts/preview_site.py"
