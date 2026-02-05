# ---------------------------------------------------------------------
# Utility functions for notebook use
# ---------------------------------------------------------------------
def get_perturbation_expression(
    mdata: MuData,
    guide_id: str,
    gene_id: str,
    guide_mod: str = "guide",
    gene_mod: str = "gene",
    gene_layer: Optional[str] = None,
    guide_assignment_layer: str = "guide_assignment",
    non_targeting_label: str = "non_targeting",
) -> pd.DataFrame:
    """
    Compare gene expression between cells with a specific guide vs non-targeting controls.

    This function identifies:
      1. "Perturbed" cells: cells assigned the specified guide
      2. "Control" cells: cells with ONLY non-targeting guides (no targeting guides)

    Then returns expression values for the specified gene in both groups.

    Parameters
    ----------
    mdata : MuData
        MuData object with guide and gene modalities.
    guide_id : str
        Guide identifier (must exist in guide.var_names).
    gene_id : str
        Gene identifier (must exist in gene.var_names).
    guide_mod : str
        Key for guide modality in mdata.mod.
    gene_mod : str
        Key for gene modality in mdata.mod.
    gene_layer : str, optional
        Layer in gene modality to use for expression. If None, uses .X.
    guide_assignment_layer : str
        Layer in guide modality containing binary assignment matrix.
    non_targeting_label : str
        Value in guide.var["label"] marking non-targeting guides.

    Returns
    -------
    pd.DataFrame
        Long-format DataFrame with columns:
          - "expression": Gene expression value
          - "group": Either the guide_id or "non_targeting"
        Index is barcode (cell ID).

    Example
    -------
    >>> expr_df = get_perturbation_expression(mdata, "CD81#strong", "ENSG00000110651")
    >>> sns.violinplot(data=expr_df, x="group", y="expression")
    """
    # Extract modalities
    guide: AnnData = mdata.mod[guide_mod]
    gene: AnnData = mdata.mod[gene_mod]

    # Validate inputs
    if guide_id not in guide.var_names:
        raise ValueError(f"guide_id '{guide_id}' not found in guide.var_names")
    if gene_id not in gene.var_names:
        raise ValueError(f"gene_id '{gene_id}' not found in gene.var_names")
    if "label" not in guide.var.columns:
        raise ValueError("guide.var must contain a 'label' column")

    # -----------------------------------------------------------------
    # Get gene expression for all cells
    # sc.get.obs_df returns DataFrame with index=cells, columns=genes
    # -----------------------------------------------------------------
    gene_expr_df = sc.get.obs_df(gene, keys=[gene_id], layer=gene_layer)
    gene_expr = gene_expr_df[gene_id]

    # -----------------------------------------------------------------
    # Identify cells with the specific guide (perturbed group)
    # guide_assignment_layer contains binary matrix: >0 means assigned
    # -----------------------------------------------------------------
    guide_assign = sc.get.obs_df(guide, keys=[guide_id], layer=guide_assignment_layer)
    cells_with_guide = guide_assign.index[guide_assign[guide_id] > 0].tolist()

    # -----------------------------------------------------------------
    # Identify cells with ONLY non-targeting guides (control group)
    # Step 1: Find cells with any non-targeting guide
    # Step 2: Find cells with any targeting guide
    # Step 3: Control = cells with NT guides but NO targeting guides
    # -----------------------------------------------------------------
    # Get non-targeting guide indices
    nt_guide_idx = guide.var.index[guide.var["label"] == non_targeting_label].tolist()
    if len(nt_guide_idx) == 0:
        raise ValueError(f"No guides with label '{non_targeting_label}' in guide.var['label']")

    # Cells with any non-targeting guide
    nt_assign = sc.get.obs_df(guide, keys=nt_guide_idx, layer=guide_assignment_layer)
    cells_with_nt = nt_assign.index[nt_assign.sum(axis=1) > 0].tolist()

    # Cells with any targeting (non-NT) guide
    targeting_guide_idx = guide.var.index[guide.var["label"] != non_targeting_label].tolist()
    targeting_assign = sc.get.obs_df(guide, keys=targeting_guide_idx, layer=guide_assignment_layer)
    cells_with_targeting = targeting_assign.index[targeting_assign.sum(axis=1) > 0].tolist()

    # Control cells: have NT guides but NO targeting guides
    cells_nt_only = list(set(cells_with_nt) - set(cells_with_targeting))

    # -----------------------------------------------------------------
    # Build output DataFrame
    # -----------------------------------------------------------------
    perturbed_expr = gene_expr.loc[cells_with_guide]
    control_expr = gene_expr.loc[cells_nt_only]

    expr_df = pd.DataFrame({
        "expression": pd.concat([perturbed_expr, control_expr]),
        "group": [guide_id] * len(perturbed_expr) + [non_targeting_label] * len(control_expr),
    })
    expr_df.index.name = "barcode"

    return expr_df