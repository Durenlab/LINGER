import numpy as np
import pandas as pd
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
    S.to_csv('cell_population_TF_RE_binding.txt',sep='\t')

## trans regulatory
## the regions
def cis_regulatory( RNA_file,ATAC_file, labels):
    import numpy as np
    import pandas as pd
    O_overlap=merge_columns_in_bed_file('Region_overlap.bed',1)
    N_overlap=merge_columns_in_bed_file('Region_overlap.bed',4)
    N_Peak=merge_columns_in_bed_file('Region.bed',1)
    O_all=merge_columns_in_bed_file('Peaks_0.bed',1)
## read the count file.
    RE=pd.read_csv('ATAC.txt',sep='\t',index_col=0)
    TG=pd.read_csv('RNA.txt',sep='\t',index_col=0)
## cell annotation
    with open('label.txt', 'r') as file:
    # Read the lines of the file and remove any trailing whitespace
        lines = [line.strip() for line in file]  
    # Create a list from the lines
        label = list(lines)
    labelset=list(set(label))
## extact the overlapped peaks.
    RE=RE.loc[N_overlap]
    from scipy.sparse import coo_matrix
    primary_s=pd.read_csv('Primary_RE_TG.txt',sep='\t')
    primary_s["RE"] = primary_s["RE"].apply(lambda x: x.split('_')[0]+':'+x.split('_')[1]+'-'+x.split('_')[2])
    primary_s = primary_s[primary_s["RE"].isin(O_overlap)]
    TGset=primary_s["TG"].unique()
    REset=O_overlap
# Create a dictionary mapping column names and row names to integer indices
    col_dict = {col: i for i, col in enumerate(TGset)}
    row_dict = {row: i for i, row in enumerate(REset)}
# Map the column names and row names to integer indices in the DataFrame
    primary_s["col_index"] = primary_s["TG"].map(col_dict)
    primary_s["row_index"] = primary_s["RE"].map(row_dict)
    # Extract the column indices, row indices, and values from the DataFrame
    col_indices = primary_s["col_index"].tolist()
    row_indices = primary_s["row_index"].tolist()
    values = primary_s["score"].tolist()
    # Create the sparse matrix using coo_matrix
    sparse_S = coo_matrix((values, (row_indices, col_indices)))
    sparse_S.colnames = TGset
    sparse_S.rownames = REset
## select the genes
    TGoverlap=list(set(TGset)&set(TG.index))
    target_col_indices = [col_dict[col] for col in TGoverlap]
    csc_S = sparse_S.tocsc()
    sparse_S = csc_S[:, target_col_indices]
    TG=TG.loc[TGoverlap]
    Dis=pd.read_csv('RE_TG_distance.txt',sep='\t',header=None)
    Dis.columns=['RE','REid','TG','TGid','dis']
    Dis["RE"] = Dis["RE"].apply(lambda x: x.split('_')[0]+':'+x.split('_')[1]+'-'+x.split('_')[2])
    Dis = Dis[Dis["RE"].isin(O_overlap)]
    Dis = Dis[Dis['TG'].isin(TGoverlap)]
    col_dict = {col: i for i, col in enumerate(TGoverlap)}
    row_dict = {row: i for i, row in enumerate(O_overlap)}
# Map the column names and row names to integer indices in the DataFrame
    Dis["col_index"] = Dis["TG"].map(col_dict)
    Dis["row_index"] = Dis["RE"].map(row_dict)
    col_indices = Dis["col_index"].tolist()
    row_indices = Dis["row_index"].tolist()
    values = Dis["dis"].tolist()
# Create the sparse matrix using coo_matrix
    sparse_dis = coo_matrix((values, (row_indices, col_indices)),shape=(len(O_overlap), len(TGoverlap)))
    sparse_dis.colnames = TGoverlap
    sparse_dis.rownames = O_overlap
    sparse_dis = sparse_dis.tocsc()
    A=sparse_dis.multiply(1 / 25000)
    A.data +=0.5
    A.data = np.exp(-A.data)
    sparse_dis=A
    TG=TG.mean(axis=1)
    TG=TG/TG.mean()+0.1
    RE=RE.mean(axis=1)
    RE=RE/RE.mean()+0.1
    sparse_S.data+=0.1
    from scipy.sparse import csc_matrix
    Score=csc_matrix(RE).T.multiply(sparse_S).multiply(sparse_dis).multiply(csc_matrix(TG))
    Score.rownames=TGoverlap
    Score.colnames=N_overlap
    coo = Score.tocoo()
# Combine the row names, column names, row indices, column indices, and values into a 3-column matrix
    combined = np.vstack((np.array(N_overlap)[coo.row], np.array(TGoverlap)[coo.col], coo.data)).T
# Save the combined matrix to a file
    np.savetxt('cell_population_cis_regulatory.txt', combined, fmt='%s',delimiter='\t')

def cell_type_specific_cis_regulatory( RNA_file,ATAC_file, labels):
    O_overlap=merge_columns_in_bed_file('Region_overlap.bed',1)
    N_overlap=merge_columns_in_bed_file('Region_overlap.bed',4)
    N_Peak=merge_columns_in_bed_file('Region.bed',1)
    O_all=merge_columns_in_bed_file('Peaks_0.bed',1)
## read the count file.
    RE=pd.read_csv('ATAC.txt',sep='\t',index_col=0)
    TG=pd.read_csv('RNA.txt',sep='\t',index_col=0)
## cell annotation
    with open('label.txt', 'r') as file:
    # Read the lines of the file and remove any trailing whitespace
        lines = [line.strip() for line in file]  
    # Create a list from the lines
        label = list(lines)
    labelset=list(set(label))
## extact the overlapped peaks.
    RE=RE.loc[N_overlap]
    from scipy.sparse import coo_matrix
    primary_s=pd.read_csv('Primary_RE_TG.txt',sep='\t')
    primary_s["RE"] = primary_s["RE"].apply(lambda x: x.split('_')[0]+':'+x.split('_')[1]+'-'+x.split('_')[2])
    primary_s = primary_s[primary_s["RE"].isin(O_overlap)]
    TGset=primary_s["TG"].unique()
    REset=O_overlap
# Create a dictionary mapping column names and row names to integer indices
    col_dict = {col: i for i, col in enumerate(TGset)}
    row_dict = {row: i for i, row in enumerate(REset)}
# Map the column names and row names to integer indices in the DataFrame
    primary_s["col_index"] = primary_s["TG"].map(col_dict)
    primary_s["row_index"] = primary_s["RE"].map(row_dict)
    # Extract the column indices, row indices, and values from the DataFrame
    col_indices = primary_s["col_index"].tolist()
    row_indices = primary_s["row_index"].tolist()
    values = primary_s["score"].tolist()
    # Create the sparse matrix using coo_matrix
    sparse_S = coo_matrix((values, (row_indices, col_indices)))
    sparse_S.colnames = TGset
    sparse_S.rownames = REset
## select the genes
    TGoverlap=list(set(TGset)&set(TG.index))
    target_col_indices = [col_dict[col] for col in TGoverlap]
    csc_S = sparse_S.tocsc()
    sparse_S = csc_S[:, target_col_indices]
    TG=TG.loc[TGoverlap]
    Dis=pd.read_csv('RE_TG_distance.txt',sep='\t',header=None)
    Dis.columns=['RE','REid','TG','TGid','dis']
    Dis["RE"] = Dis["RE"].apply(lambda x: x.split('_')[0]+':'+x.split('_')[1]+'-'+x.split('_')[2])
    Dis = Dis[Dis["RE"].isin(O_overlap)]
    Dis = Dis[Dis['TG'].isin(TGoverlap)]
    col_dict = {col: i for i, col in enumerate(TGoverlap)}
    row_dict = {row: i for i, row in enumerate(O_overlap)}
# Map the column names and row names to integer indices in the DataFrame
    Dis["col_index"] = Dis["TG"].map(col_dict)
    Dis["row_index"] = Dis["RE"].map(row_dict)
    col_indices = Dis["col_index"].tolist()
    row_indices = Dis["row_index"].tolist()
    values = Dis["dis"].tolist()
# Create the sparse matrix using coo_matrix
    sparse_dis = coo_matrix((values, (row_indices, col_indices)),shape=(len(O_overlap), len(TGoverlap)))
    sparse_dis.colnames = TGoverlap
    sparse_dis.rownames = O_overlap
    sparse_dis = sparse_dis.tocsc()
    A=sparse_dis.multiply(1 / 25000)
    A.data +=0.5
    A.data = np.exp(-A.data)
    sparse_dis=A
    from scipy.sparse import csc_matrix
    for label0 in labelset:
        TG_temp=TG.values[:,np.array(label)==label0].mean(axis=1)
        TG_temp=TG_temp/TG_temp.mean()+0.1
        RE_temp=RE.values[:,np.array(label)==label0].mean(axis=1)
        RE_temp=RE_temp/RE_temp.mean()+0.1
        sparse_S.data+=0.1
        Score=csc_matrix(RE_temp).T.multiply(sparse_S).multiply(sparse_dis).multiply(csc_matrix(TG_temp))
        coo = Score.tocoo()
# Combine the row names, column names, row indices, column indices, and values into a 3-column matrix
        combined = np.vstack((np.array(N_overlap)[coo.row], np.array(TGoverlap)[coo.col], coo.data)).T
# Save the combined matrix to a file
        np.savetxt('cell_type_specific_cis_regulatory_'+label0+'.txt', combined, fmt='%s',delimiter='\t')