import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import pandas as pd
import numpy as np
import seaborn as sns
import plotly.express as px
from plotly.subplots import make_subplots

def venn_diagram(label1, label1_size, label2, label2_size, both_size, size_label='', overlap=False):
    """
    Generate a Venn diagram to visualize the relationships between two sets.

    Parameters:
    label1 (str): Label for the first set.
    label1_size (int): Size of the first set.
    label2 (str): Label for the second set.
    label2_size (int): Size of the second set.
    both_size (int): Size of the overlap between the two sets.
    size_label (str, optional): Label to indicate the unit of sizes, e.g., 'MB', 'GB'.
    overlap (bool, optional): If True, adjust set sizes to account for the overlap region.

    Returns:
    the figure

    This function generates a Venn diagram using the matplotlib_venn library to visualize the relationships
    between two sets and their overlap. The function allows for adjusting set sizes when there is an overlap,
    by subtracting the overlap size from the respective set sizes before plotting.
    """

    if overlap:
        label1_size = label1_size - both_size
        label2_size = label2_size - both_size

    # Create the Venn diagram
    plt.figure(figsize=(16, 10))
    font1 = {'size': 15, 'style': 'italic'}
    plt.rc('font', **font1)
    p = venn2(subsets=(label1_size, label2_size, both_size),
              set_labels=(f'{label1} \n {label1_size + both_size:,}{size_label}',
                          f'{label2} \n {label2_size + both_size:,}{size_label}'))

    # Update the text for set labels
    for text in p.set_labels:
        print(text.get_text())
        # You can modify the font size and other properties here if needed
        # text.set_fontsize(18)
        # text.set_text(text.get_text())

    # Update the text for overlap region labels
    for i in ['10', '11', '01']:
        p.get_label_by_id(i).set_text(f"{int(p.get_label_by_id(i).get_text()):,}{size_label}")

    return p



CLASS_IDX_START = 10

def plot_distribution_with_sd(data_dict, title, x_label):
    """
    Create a distribution plot with signs for mean and standard deviation using Seaborn.
    The plot is cut at the 98th percentile and the area above 95th percentile is colored.
    
    Parameters:
    data_dict (dict): A dictionary containing {'label': data} pairs, where data is an array-like of numerical data.
    title (str): Title for the plot.
    """
    plt.figure(figsize=(8, 6))
    
    # Calculate the 98th percentile to cut the plot
    all_data = np.concatenate(list(data_dict.values()))
    percentile_98 = np.percentile(all_data, 98)
    percentile_95 = np.percentile(all_data, 95)
    
    colors = sns.color_palette('pastel')
    colors2 = sns.color_palette('Set1')
    
    for idx, (label, data) in enumerate(data_dict.items()):
        sns.histplot(data, kde=True, label=label, color=colors[idx], alpha=0.6, stat ='density')
    
        mean_value = np.mean(data)
        std_value = np.std(data)
    
        plt.axvline(mean_value, color=colors2[len(colors2) - (idx + 1) ], linewidth=1, label=f'{label} Mean')
        plt.axvline(mean_value + std_value, color=colors[idx], linestyle='dashed', linewidth=1, label=f'{label} 1 SD')
        plt.axvline(mean_value - std_value, color=colors[idx], linestyle='dashed', linewidth=1)
    
    plt.xlim(0, percentile_98)  # Set x-axis limit up to the 98th percentile
    plt.xlabel(x_label)
    plt.ylabel('Frequency')
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.show()

    

def plot_scatter_with_info_and_histogram(df, x_data, y_data, info_data,title):
    """
    Create a scatter plot with additional information and a histogram using Plotly.
    
    Parameters:
    df (DataFrame): The DataFrame containing the data.
    x_data (str): Column name for the x-axis data.
    y_data (str): Column name for the y-axis data.
    info_data (str): Column name for additional information to display on the plot.
    """
    fig = make_subplots(rows=2, cols=1, column_widths=[1], row_heights=[0.7,0.3], shared_yaxes=True,
                        shared_xaxes=True,
                        horizontal_spacing=0.05)#, subplot_titles=("Scatter Plot", "Histogram"))
    
    scatter_fig = px.scatter(df, x=x_data, y=y_data, hover_name=info_data)
    scatter_fig.update_traces(textposition='top center')
    scatter_fig.update_xaxes(title_text=x_data, title_font=dict(size=14), tickangle=45)
    scatter_fig.update_yaxes(title_text=y_data, title_font=dict(size=14))
    

    
    histogram_fig = px.histogram(df, x=x_data)
    histogram_fig.update_xaxes(title_text="Variants Count", title_font=dict(size=14))
    histogram_fig.update_yaxes(title_text="Frequency", title_font=dict(size=14))
    
    scatter_trace = scatter_fig['data'][0]
    histogram_trace = histogram_fig['data'][0]
    
    fig.add_trace(scatter_trace, row=1, col=1)
    fig.add_trace(histogram_trace, row=2, col=1)
    
    fig.update_layout(title_text=title,
                      title_font=dict(size=16),
                      showlegend=False)
    
    fig.show()

    
def get_scatter_plot(class_score, diff_threshold, title="" ):
   
    name = class_score.chrom + '_' + class_score.pos.astype(str) \
    + class_score.ref + class_score.alt
    name = name.rename('name')
    max_diff_class =np.abs(class_score.iloc[:,CLASS_IDX_START:]).idxmax(axis=1).rename('class')
    scatter_plot_df = pd.concat([name,class_score.seqclass_max_absdiff, max_diff_class], axis=1)
    plot_scatter_with_info_and_histogram(scatter_plot_df,'class','seqclass_max_absdiff', 'name', title )
    

def described_variant_hist(name_col,title='# of described variants'):
    name_col_bool = name_col == '.'
    mask = name_col_bool.copy()
    name_col_bool.loc[mask] = 'novel'
    name_col_bool.loc[~mask] = 'known'
    px.histogram(name_col_bool, title=title)
    return px.histogram(name_col_bool, title=title)


def plot_box_plot(series,title,x_label, without_outliers=True):
    # Calculate the quartiles for the series
    q1 = series.quantile(0.25)
    q3 = series.quantile(0.75)
    iqr = q3 - q1

    if without_outliers:
        # Define the lower and upper bounds for non-outliers
        lower_bound = q1 - 1.5 * iqr
        upper_bound = q3 + 1.5 * iqr

        # Filter out outliers
        series = series[(series >= lower_bound) & (series <= upper_bound)]

    # Create a Plotly box plot without outliers
    fig = px.box(series, x=series)

    # Update the x-axis label
    fig.update_layout(xaxis_title=x_label)
    fig.update_layout(title=title)


    # Show the plot
    fig.show()