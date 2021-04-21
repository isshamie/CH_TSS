import plotly.graph_objects as go
import pandas as pd
import numpy as np
from mplh.color_utils import get_colors
#from tss.visualize.co import


def create_sankey(df, ordered_cols, init="TSS", source=None, out_df=(),
                  is_init=True, add_colname=False,
                  var_name=None):  # , node_color_dict=None):
    if source is None:  # initial
        out_df = []
        source = init
    elif len(ordered_cols) == 0:  # Stop criteria
        return out_df  # pd.DataFrame(out_df, columns=["source", "target", "value"])
    curr_col = ordered_cols[0]
    for target, val in df.groupby(curr_col):
        if var_name is None:
            var = val.shape[0]
        else:
            var = val[var_name].sum()

        if add_colname:
            out_df.append((source, curr_col + f"_{target}", var))
            new_source = curr_col + f"_{target}"
        else:
            out_df.append((source, target, var))
            new_source = target
        out_df = create_sankey(val, ordered_cols[1:], source=new_source,
                               out_df=out_df, is_init=False,
                               add_colname=add_colname,
                               var_name=var_name)
    if is_init:
        return pd.DataFrame(out_df,
                            columns=["source", "target", "value"])
    return out_df


def array_to_rgb(x, has_opac=False):
    curr = x[:3]
    curr = curr.astype(int).astype(str)
    out = "rgb(%s)" % (",".join(curr))
    if has_opac:
        out = out.replace(")", ",%0.1f)" % x[-1])
        out = out.replace('rgb', 'rgba')
    return out


def plot_sankey(link_dat, label_dict, label_colors, out_f=None):
    fig = go.Figure(data=[go.Sankey(# Define nodes
        node=dict(pad=15, thickness=15,
            line=dict(color="black", width=0.5),
            label=list(label_dict.keys()), color=label_colors),
        link=link_dat, )])

    fig.update_layout(title_text="TSSs Detected", font_size=10)
    #fig.show()

    if out_f is not None:
        fig.write_image(out_f + ".pdf")
        fig.write_image(out_f + ".png")
    return


def wrap_sankey(df, ordered_cols, var_name=None, node_opacity=0.8,
                link_opacity=0.4, name="TSSs detected",
                add_colname=True, out_f=None):

    sankey_df = create_sankey(df, ordered_cols=ordered_cols,
                              add_colname=add_colname, init=name,
                              var_name=var_name)
    label_dict = {}
    for i in ['source', 'target']:
        for j in sankey_df[i]:
            if j not in label_dict:
                label_dict[j] = len(label_dict)
    print(label_dict)
    sankey_df = sankey_df.replace(label_dict, )

    # Process the colors
    # For each source, get the np array colors, convert to rgba, add opacity, and make the links lower opacity
    list(label_dict.keys())
    label_colors = get_colors(scheme='categorical',
                              n_colors=len(label_dict), use_white=False,
                              use_black=False, return_p=True)
    label_colors = label_colors[2] * 255
    label_colors = np.floor(label_colors).astype(int)
    label_colors = np.append(label_colors, node_opacity * np.ones(
        [label_colors.shape[0], 1]), axis=1)
    label_colors_link = label_colors.copy()  # [np.append(x, opacity) for x in label_colors]
    label_colors_link = np.append(label_colors, link_opacity * np.ones(
        [label_colors.shape[0], 1]), axis=1)
    label_colors = ([array_to_rgb(x, True) for x in label_colors])
    label_colors_link = (
    ([array_to_rgb(x, True) for x in label_colors_link]))

    # Add opacity for the links to be less
    sankey_df["color"] = sankey_df["source"].apply(
        lambda x: label_colors_link[x])

    link_dat = {x: list(sankey_df[x].values) for x in
                sankey_df.to_dict()}
    # Plot
    plot_sankey(link_dat, label_dict, label_colors, out_f=out_f)

    # put names back in place
    sankey_df[["source", "target"]] = sankey_df[
        ['source', 'target']].replace(
        {label_dict[x]: x for x in label_dict})
    if out_f is not None:
        # Save csv file
        sankey_df.sort_values("value").to_csv(out_f+".csv")
    return sankey_df, label_colors, label_dict


def get_venn(peaks, group_A_list, group_B_list, return_sets=False):
    group_A_peaks = set()
    for i in group_A_list:
        group_A_peaks = group_A_peaks.union(peaks[i])

    group_B_peaks = set()
    for i in group_B_list:
        group_B_peaks = group_B_peaks.union(peaks[i])

    if return_sets:
        return group_A_peaks, group_B_peaks
    else:
        # Venn Diagram
        a = group_A_peaks - group_B_peaks
        b = group_B_peaks - group_A_peaks
        overlap = group_B_peaks.intersection(group_A_peaks)
        return a, b, overlap

