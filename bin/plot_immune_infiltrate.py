#!/usr/bin/env python3
import argparse
import pandas as pd
import matplotlib.pyplot as plt


def get_data_to_plot(immune_panel_results_df):
    # Filter rows based on column: 'Biomarker Type'
    immune_panel_results_df = immune_panel_results_df[immune_panel_results_df['Biomarker Type'] == "immune_inf"]
    # Derive column 'Cell Name' from column: 'Biomarker name'
    immune_panel_results_df["Cell Name"] = immune_panel_results_df["Biomarker name"].str.split("_").str[0]
    # Select columns: 'Cell Name', 'Result'
    immune_panel_results_df = immune_panel_results_df.loc[:, ['Cell Name', 'Result']]
    # Sort by column: 'Result' (ascending)
    immune_panel_results_df = immune_panel_results_df.sort_values(['Result'], ascending=False)
    return immune_panel_results_df


def main(args):
    
    celltype_colors = {
    "Cancer": "#1f77b4",
    "Fibroblast": "#7f7f7f",
    "NK": "#2ca02c",
    "CD4": "#d62728",
    "Treg": "#9467bd",
    "Endothelial": "#8c564b",
    "Bcell": "#e377c2",
    "Eosinophil": "#ff7f0e",
    "Neutrophil": "#bcbd22",
    "CD8": "#17becf",
}
    
    immune_panel_results_df = pd.read_csv(args.result_data)
    
    plot_data_df = get_data_to_plot(immune_panel_results_df.copy())
    
    # Get cell type of any Nan vals
    missing_cell_types = plot_data_df[plot_data_df['Result'].isna()]['Cell Name'].tolist()
    # remove nan values from plot_data_df
    plot_data_df = plot_data_df[plot_data_df['Result'].notna()]
    
    # Change names to include percentage
    plot_data_df.insert(1, "Cell Type", plot_data_df.apply(lambda row : row["Cell Name"] + " (" + f'{row["Result"] * 100:01.0f}%' + ")", axis=1))
    plot_data_df.drop(columns=['Cell Name'], inplace=True)
    
    # Transpose data for plotting
    plot_data_df_T = plot_data_df.T
    plot_data_df_T.columns = plot_data_df_T.iloc[0]
    plot_data_df_T = plot_data_df_T[1:]
    
    # get colours
    colours_list = [celltype_colors[c.split(' ')[0]] for c in plot_data_df_T.columns]
    ax = plot_data_df_T.plot.bar(stacked=True, rot=0, figsize=(2, 6), color=colours_list)

    # legend order reversed and place legend outside plot
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], bbox_to_anchor=(1.04, 0.75), loc='upper left', borderaxespad=0)
    # plt.subplots_adjust(right=0.25)  # increase/decrease 0.75 to fit legend
    
    plt.savefig(args.out, dpi=300, bbox_inches='tight', pad_inches=0.1)


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--result_data', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    
    main(args)
