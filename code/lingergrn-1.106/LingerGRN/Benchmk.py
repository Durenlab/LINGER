import seaborn as sns
import matplotlib.colors as mcolors
def generate_colors(N):
    """
    Generate N visually appealing colors using seaborn color palette.
    
    Args:
        N (int): The number of colors to generate.
    
    Returns:
        list: A list of N RGB tuples representing the generated colors.
    """
    color_palette = sns.color_palette("husl", N)
    colors = [mcolors.rgb2hex(color_palette[i]) for i in range(N)]
    return colors

def bm_trans(TFName,Method_name,Groundtruth,Infer_trans,outdir,filetype):
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve, roc_auc_score
    import seaborn as sns
# Set the working directory
# Use the appropriate path for your system
    import os
    data0=pd.read_csv(Groundtruth,sep='\t', skiprows=5,header=0)
    data1=data0.groupby(['symbol'])['score'].max()
    label=data1.sort_values(axis=0, ascending=False)
    label=label.reset_index()
    N = 1000  # top 500 TG as the ground truth
    label = label['symbol'].iloc[:N]
    # Load data
    colors=generate_colors(len(Infer_trans))
    for i in range(len(Infer_trans)):
        if filetype=='list':
            data2 = pd.read_csv(Infer_trans[i], sep='\t')
        #label = pd.read_csv('~/SC_NET/all_data/result/LL_sc/max/' + info0['id'][i] + '_gene_score_5fold.txt', header=True)
            loc = data2[data2['TF'] == TFName].index
            TGName = data2['TG'].values[loc].tolist()
            TGset = TGName.copy()
            Score = data2['score'].values[loc]
        if filetype=='matrix':
            data2=pd.read_csv(Infer_trans[i], sep='\t',header=0,index_col=0)
            TGset=data2.index
            Score=data2[TFName].values
        d1 = np.zeros(len(TGset))
        loc = np.where(np.isin(TGset, label))[0]
        d1[loc] = 1
        fpr, tpr, thresholds = roc_curve(d1, Score)
# Compute AUC
        auc = roc_auc_score(d1, Score)
    # Plot the ROC curve
        plt.plot(fpr, tpr, color=colors[i], label=Method_name[i]+': (AUC = %0.2f)' % auc)
        plt.plot([0, 1], [0, 1], color='black', linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('Receiver Operating Characteristic')
        plt.legend(loc="lower right")
# Display the plot in the notebook
    plt.savefig(outdir+"trans_roc_curve"+TFName+".png", format='png', bbox_inches='tight')
    plt.show()
    plt.close()
    from sklearn.metrics import precision_recall_curve
    from sklearn.metrics import average_precision_score
    for i in range(len(Infer_trans)):
        if filetype=='list':
            data2 = pd.read_csv(Infer_trans[i], sep='\t')
        #label = pd.read_csv('~/SC_NET/all_data/result/LL_sc/max/' + info0['id'][i] + '_gene_score_5fold.txt', header=True)
            loc = data2[data2['TF'] == TFName].index
            TGName = data2['TG'].values[loc].tolist()
            TGset = TGName.copy()
            Score = data2['score'].values[loc]
        if filetype=='matrix':
            data2=pd.read_csv(Infer_trans[i], sep='\t',header=0,index_col=0)
            TGset=data2.index
            Score=data2[TFName].values
        d1 = np.zeros(len(TGset))
        loc = np.where(np.isin(TGset, label))[0]
        d1[loc] = 1
# Assuming you have the true labels (y_true) and predicted probabilities (y_scores) for your classifier
# Calculate the average precision score (AUPR
        aupr = average_precision_score(d1, Score)
        auprr=aupr*len(d1)/sum(d1)
# Assuming you have the true labels (y_true) and predicted probabilities (y_scores) for your classifier
# Calculate precision and recall values
        precision, recall, _ = precision_recall_curve(d1, Score)
    # Plot precision-recall curve
        plt.plot(recall, precision, color=colors[i], label=Method_name[i]+': (AUPR ratio = %0.2f)' % auprr)
        plt.ylim([0.0, 0.6])
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.title('Precision-Recall Curve')
        plt.legend(loc='upper right')
    plt.savefig(outdir+"trans_pr_curve"+TFName+".png", format='png', bbox_inches='tight')
    plt.show()
    plt.close()