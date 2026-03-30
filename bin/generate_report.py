#!/usr/bin/env python3
import argparse
import pandas as pd
from jinja2 import Template

parser = argparse.ArgumentParser()
parser.add_argument('--panel', required=True)
parser.add_argument('--template', required=True)
parser.add_argument('--scores', required=True)
parser.add_argument('--out', required=True)
parser.add_argument('--sample', required=True)
args = parser.parse_args()

panel_df = pd.read_csv(args.panel, index_col=0, dtype={'ID': str})
scores_df = pd.read_csv(args.scores, index_col=0)

all_results_df = scores_df.merge(panel_df[['Biomarker name','Biomarker Type','Result Options']], left_index=True, right_index=True, how='left')
all_results_df.reset_index(inplace=True)

with open(args.template, 'r') as fh:
    template = Template(fh.read(), trim_blocks=True)

results = {
    'sample_name': args.sample,
    'panel_name': '',
    'url': 'https://github.com/lucy924/nanopore_multiBM_pipeline',
    'all_results': all_results_df.values.tolist()
}

rendered = template.render(repo=results)
with open(args.out, 'w') as fw:
    fw.write(rendered)
