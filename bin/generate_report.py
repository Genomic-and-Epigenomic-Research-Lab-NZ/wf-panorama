#!/usr/bin/env python3

import argparse
from jinja2 import Template
import codecs
import pandas as pd
import markdown
from shared_functions import BIOMARKER_ID, BIOMARKER_NAME, get_BM_TYPE_FULL, BIOMARKER_TYPE, RESULT_OPTIONS, RESULT


def main(args):
    
    BIOMARKER_TYPE_FULL = get_BM_TYPE_FULL(path2panel=args.panel)

    # log = open(args.log, 'w')
    
    # create a dict with all data that will populate the template
        
    # Get results from csv into dict
    results_df = pd.read_csv(args.scores, index_col=0)
    panel_df = pd.read_csv(args.panel, index_col=0, dtype={BIOMARKER_ID: str})
    panel_df.rename(columns = {BIOMARKER_TYPE_FULL: BIOMARKER_TYPE}, inplace=True)
        
    header = [BIOMARKER_ID, BIOMARKER_NAME, BIOMARKER_TYPE, RESULT_OPTIONS, RESULT, "Score"]
    panel_cols = [BIOMARKER_NAME, BIOMARKER_TYPE, RESULT_OPTIONS]

    # Merge using index
    all_results_df = results_df.merge(panel_df[panel_cols], left_index=True, right_index=True, how='left', )
    all_results_df.reset_index(names='ID', inplace=True)
    all_results_df = all_results_df[header]
        
    # all_results_dict = all_results_df.to_dict()
    all_results_df = all_results_df.map(lambda x: x.replace('|', '\\|') if isinstance(x, str) else x)

    results = {
        "report_title": args.report_title,
        "sample_name": args.sample,
        "panel_name": args.project_name,
        "url": "https://github.com/Genomic-and-Epigenomic-Research-Lab-NZ/wf-panorama"
    }

    all_data_lists = list()
    for idx, entry in all_results_df.iterrows():
        # entry_list = [idx]
        # entry_list.extend(list(entry))
        all_data_lists.append(list(entry))
        
    results["all_results"] = all_data_lists

    # get word result from normalised score
    normalised_score = all_results_df['Score'].iloc[-1]
    if normalised_score < 25:
        results["word_result"] = 'Highly Likely'
    elif normalised_score < 50:
        results["word_result"] = 'Somewhat Likely'
    elif normalised_score < 75:
        results["word_result"] = 'Somewhat Unlikely'
    else:
        results["word_result"] = 'Highly Unlikely'
        

    # render the template

    with open(args.template, "r") as file:
        template = Template(file.read(), trim_blocks=True)
    rendered_file = template.render(repo=results)

    # output the markdown file
    output_file = codecs.open(args.out, "w", "utf-8")
    output_file.write(rendered_file)
    output_file.close()

    # output the html file
    html_body = markdown.markdown(rendered_file, extensions=['tables'])
    html_content = f"""<!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>{args.report_title}</title>
        </head>
        <body>
        {html_body}
        </body>
        </html>"""
    
    html_out = args.out.rsplit('.', 1)[0] + '.html'
    output_html = codecs.open(html_out, "w", "utf-8")
    output_html.write(html_content)
    output_html.close()

    # log.close()


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--panel', required=True)
    parser.add_argument('--template', required=True)
    parser.add_argument('--scores', required=True)
    parser.add_argument('--report_title', required=True, default='Report Title')
    parser.add_argument('--project_name', required=True)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--log', required=True)
    args = parser.parse_args()
    
    main(args)

# Other weird stuff
# panel_df = pd.read_csv(args.panel, index_col=0, dtype={'ID': str})
# scores_df = pd.read_csv(args.scores, index_col=0)

# all_results_df = scores_df.merge(panel_df[['Biomarker name','Biomarker Type','Result Options']], left_index=True, right_index=True, how='left')
# all_results_df.reset_index(inplace=True)

# with open(args.template, 'r') as fh:
#     template = Template(fh.read(), trim_blocks=True)

# results = {
#     'sample_name': args.sample,
#     'panel_name': '',
#     'url': 'https://github.com/lucy924/nanopore_multiBM_pipeline',
#     'all_results': all_results_df.values.tolist()
# }

# rendered = template.render(repo=results)
# with open(args.out, 'w') as fw:
#     fw.write(rendered)
