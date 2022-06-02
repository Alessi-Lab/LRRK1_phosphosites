import numpy as np
import pandas as pd
import re
import seaborn as sns
import matplotlib.pylab as plt


uniprot = re.compile("(?P<accession>[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(?P<isotype>-\d+)?")

proteases = ["AspN",
                 #"Chymotrypsin",
                 "ChyTryp",
                 "Tryp"
                 ]

if __name__ == "__main__":
    results = {}
    input_file = r"C:\Users\toanp\Downloads\LRRK1_MQ (1).txt"
    intensity_df = pd.read_csv(input_file, sep="\t")
    fas = r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\For Phospho Motif\Updated\Human_UniprotSP_Cano+Iso_052021.fasta"
    seqs = {}

    with open(fas, "rt") as fasFile:
        current_seq = ""
        current_id = ""
        for line in fasFile:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    if current_seq:
                        seqs[current_id[:]] = current_seq[:]

                acc = uniprot.search(line)
                if acc:
                    current_seq = ""
                    current_id = acc.group(0)
            else:
                current_seq += line
    merged_df = []
    for i, r in intensity_df.iterrows():
        pos = int(r["Positions within proteins"])-1
        seq = \
            seqs[r["Proteins"]][pos-6:pos] + \
            seqs[r["Proteins"]][pos].lower() + \
            seqs[r["Proteins"]][pos+1:pos+6]
        for c in intensity_df.columns[0:40]:
            if c.startswith("Intensity "):
                for p in proteases:
                    if p in c:
                        merged_df.append([
                            pos+1,
                            p,
                            c.replace("Intensity ", ""),
                            r["Localization prob"],
                            seq,
                            r[c]
                        ])
                        break
    merged_df = pd.DataFrame(merged_df, columns=["Position", "Proteases", "Sample", "Score", "Window", "Abundance"])
    merged_df.sort_values("Position", inplace=True)
    e = 1
    n = 500

    samples = [c.replace("Intensity ", "") for c in intensity_df.columns[28:40]]
    samples_columns = []
    condition_dict = {}
    for p in proteases:
        if p not in condition_dict:
            condition_dict[p] = []
        for s in samples:
            if p in s.split("_"):
                condition_dict[p].append(s)

    for p in proteases:
        for i in range(len(condition_dict[p]), 0, -2):
            for s in condition_dict[p][i-2:i]:
                samples_columns.append((p, s))
        #for s in condition_dict[p]:
            #samples_columns.append((p, s))


    multiindex = pd.MultiIndex.from_tuples(samples_columns, names=["Proteases", "Sample"])

    c = merged_df
    fontsize_pt = plt.rcParams['ytick.labelsize']
    dpi = 72.27
    top_margin = 0.2
    bottom_margin = 0.2
    left_margin = 0.2
    right_margin = 0.2
    figure_height = (len(c.index)/30) / (1 - top_margin - bottom_margin)
    figure_width = 10 / (1-left_margin-right_margin)
    c = c.set_index([
        #"Phospho",
        "Position", "Sample", "Proteases", "Window"])
    c = c.unstack("Proteases")

    b = pd.pivot_table(c, values="Abundance", columns="Sample", index=["Position",
                                                                        #"Phospho",
                                                                        "Window"])
    b.fillna(0, inplace=True)
    b = b.T
    for i in b.columns:
        b0 = b[i][b[i]==0]
        b[i] = (np.log2(b[i], where=b[i]>0) - np.log2(b[i], where=b[i]>0).mean()) / np.log2(b[i], where=b[i]>0).std(ddof=1)
        for ind in b0.index:
            b[i].loc[ind] = np.nan
    b = b.T
    new_df = pd.DataFrame(index=b.index, columns=multiindex)
    for i in new_df.columns:
        if i in b.columns:
            new_df[i] = b[i]
        else:
            new_df[i].fillna(0, inplace=True)

    new_df.to_csv(input_file + f".reversed.csv")
    fig, ax = plt.subplots(
        figsize=(figure_width, figure_height),
        gridspec_kw=dict(top=1-top_margin, bottom=bottom_margin, left=left_margin, right=1-right_margin)
    )
    mask = np.isnan(b)

    sns.heatmap(new_df, cmap="YlGnBu", square=True, ax=ax)
    ax.set_facecolor("silver")
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    for label in ax.get_yticklabels():
        label.set_weight("bold")
    for label in ax.get_xticklabels():
        label.set_weight("bold")
    plt.xticks(rotation=90)
    plt.savefig(input_file + f"result.reversed.pdf")

