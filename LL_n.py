def merge_columns_in_bed_file(file_path,startcol):
    merged_values = []
    
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            col1 = columns[-1+startcol]
            col2 = columns[startcol]
            col3 = columns[1+startcol]
            merged_value = f"{col1}:{col2}-{col3}"
            merged_values.append(merged_value)
    return merged_values
def merge_columns_in_bed_file2(file_path,startcol):
    merged_values = []
    
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            col1 = columns[-1+startcol]
            col2 = columns[startcol]
            col3 = columns[1+startcol]
            merged_value = f"{col1}_{col2}_{col3}"
            merged_values.append(merged_value)
    return merged_values
def format_RE_tran12(region):
    chr, range_ = region.split(":")
    start, end = range_.split("-") 
    return "_".join([chr, start, end])

def cell_type_specific_TF_RE_binding( RNA_file,ATAC_file, labels):
    ## the regions
    O_overlap=merge_columns_in_bed_file('Region_overlap.bed',1)
    N_overlap=merge_columns_in_bed_file('Region_overlap.bed',4)
    N_Peak=merge_columns_in_bed_file('Region.bed',1)
    O_all=merge_columns_in_bed_file('Peaks_0.bed',1)
    import numpy as np
    import pandas as pd
## read the count file.
    RE=pd.read_csv(ATAC_file,sep='\t',index_col=0)
    TG=pd.read_csv(RNA_file,sep='\t',index_col=0)
## cell annotation
    with open(labels, 'r') as file:
    # Read the lines of the file and remove any trailing whitespace
        lines = [line.strip() for line in file]  
    # Create a list from the lines
        label = list(lines)
    labelset=list(set(label))

## extact the overlapped peaks.
    RE=RE.loc[N_overlap]
    TFbinding=pd.read_csv('TF_binding_f.txt',sep='\t',index_col=0)
    TFbinding = TFbinding.loc[O_overlap]
    mat=pd.read_csv('Primary_TF_RE.txt',sep='\t',index_col=0)
    TFs = mat.columns
    TFoverlap = list(set(TFs) & set(TG.index))
    mat = mat.loc[O_overlap] 
    mat=mat[TFoverlap]
    TFbinding = TFbinding[TFoverlap]
    TF = TG.loc[TFoverlap]
    mat_m=np.mean(mat.values[mat>0])
    mat = mat / mat_m
    mat.values[mat.values<0]=0
    TFbinding = TFbinding / TFbinding.mean(axis=1).mean()
    labelset=list(set(label))
    for label0 in labelset:
        TF_cluster = TF.values[:,np.array(label)==label0].mean(axis=1)
        TF_cluster = TF_cluster[None,:]
        RE_cluster = RE.values[:,np.array(label)==label0].mean(axis=1)
        RE_cluster = RE_cluster[:,None]
        S = np.log(RE_cluster+0.1) + np.log(mat+TFbinding+0.1) + np.log(TF_cluster+0.1)
        S.index=N_overlap
        S.to_csv('cell_type_specific_TF_RE_binding'+label0+'.txt',sep='\t')

def TF_RE_binding( RNA_file,ATAC_file, labels):
    ## the regions
    O_overlap=merge_columns_in_bed_file('Region_overlap.bed',1)
    N_overlap=merge_columns_in_bed_file('Region_overlap.bed',4)
    N_Peak=merge_columns_in_bed_file('Region.bed',1)
    O_all=merge_columns_in_bed_file('Peaks_0.bed',1)
    import numpy as np
    import pandas as pd
## read the count file.
    RE=pd.read_csv(ATAC_file,sep='\t',index_col=0)
    TG=pd.read_csv(RNA_file,sep='\t',index_col=0)
## cell annotation
    with open(labels, 'r') as file:
    # Read the lines of the file and remove any trailing whitespace
        lines = [line.strip() for line in file]  
    # Create a list from the lines
        label = list(lines)
    labelset=list(set(label))
## extact the overlapped peaks.
    RE=RE.loc[N_overlap]
    TFbinding=pd.read_csv('TF_binding_f.txt',sep='\t',index_col=0)
    TFbinding = TFbinding.loc[O_overlap]
    mat=pd.read_csv('Primary_TF_RE.txt',sep='\t',index_col=0)
    TFs = mat.columns
    TFoverlap = list(set(TFs) & set(TG.index))
    mat = mat.loc[O_overlap] 
    mat=mat[TFoverlap]
    TFbinding = TFbinding[TFoverlap]
    TF = TG.loc[TFoverlap]
    mat_m=np.mean(mat.values[mat>0])
    mat = mat / mat_m
    mat.values[mat.values<0]=0
    TFbinding = TFbinding / TFbinding.mean(axis=1).mean()
    labelset=list(set(label))
    TF_cluster = TF.values.mean(axis=1)
    TF_cluster = TF_cluster[None,:]
    RE_cluster = RE.values.mean(axis=1)
    RE_cluster = RE_cluster[:,None]
    S = np.log(RE_cluster+0.1) + np.log(mat+TFbinding+0.1) + np.log(TF_cluster+0.1)
    S.index=N_overlap
    S.to_csv('cell_population_TF_RE_binding_'+label0+'.txt',sep='\t')