import argparse
import pandas as pd
import numpy as np
from anndata import AnnData
from sklearn.preprocessing import RobustScaler
from pycytominer import feature_select


def scale_per_group(
    adata,
    groupby_keys=['plate', 'percent_volume_dmso'],
    control_label='DMSO',
    source_layer=None,
    target_layer='scaled',
    strict=False
):
    X = adata.X if source_layer is None else adata.layers[source_layer]
    X = X.copy()
    obs = adata.obs.copy()
    scaled = np.full(X.shape, np.nan, dtype=np.float32)
    group_labels = obs[groupby_keys].astype(str).agg('|'.join, axis=1)

    for group in group_labels.unique():
        group_mask = group_labels == group
        control_mask = group_mask & (obs['compound'] == control_label)

        if np.sum(control_mask) < 2:
            if strict:
                raise ValueError(f"Insufficient controls in group: {group}")
            continue

        scaler = RobustScaler().fit(X[control_mask])
        scaled[group_mask] = scaler.transform(X[group_mask])

    adata.layers[target_layer] = scaled
    return adata


def apply_feature_selection(
    adata,
    layer='scaled',
    operations=['variance_threshold', 'correlation_threshold'],
    copy=True
):
    if copy:
        adata = adata.copy()

    X = adata.to_df(layer=layer)
    feature_columns = X.columns.tolist()
    X_selected = feature_select(
        X,
        features=feature_columns,
        operation=operations
    )
    adata = adata[:, X_selected.columns]
    return adata


def main(metadata_file, features_file, output_file, groupby_keys, control_label):
    print("📥 Loading metadata and feature tables...")
    obs = pd.read_csv(metadata_file).set_index('sample_id')
    X = pd.read_csv(features_file).drop('container_id', axis=1).set_index('sample_id')

    print(f"🔎 Metadata shape: {obs.shape}")
    print(f"🔎 Features shape: {X.shape}")

    idx = obs.index.intersection(X.index)
    obs = obs.reindex(idx)
    X = X.reindex(idx)
    obs.index = obs.index.astype(str)  # 👈 fix warning
    print(f"🔗 Aligned data on {len(idx):,} sample_ids")

    print("📦 Constructing AnnData object...")
    adata = AnnData(X=X, obs=obs)
    print(f"🧪 AnnData shape: {adata.shape}")

    print("⚖️ Scaling features per group using RobustScaler...")
    adata = scale_per_group(
        adata,
        groupby_keys=groupby_keys,
        control_label=control_label
    )
    print(f"✅ Scaling complete. Scaled data stored in `adata.layers['scaled']`.")

    print("🧹 Running feature selection via pycytominer...")
    adata = apply_feature_selection(adata, copy=False)
    print(f"✅ Feature selection complete. Final shape: {adata.shape}")

    print("🧼 Cleaning obs metadata to ensure h5ad compatibility...")
    adata.obs = adata.obs.astype(str)
    print("✅ Metadata cleaned.")

    print(f"💾 Writing processed AnnData to {output_file}...")
    adata.write(output_file)
    print("🎉 Done.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process Cell Painting data")
    parser.add_argument("--metadata", required=True, help="Path to metadata CSV file")
    parser.add_argument("--features", required=True, help="Path to features CSV file")
    parser.add_argument("--output", required=True, help="Output .h5ad file path")
    parser.add_argument(
        "--groupby",
        nargs="+",
        default=["container_id", "percent_volume_dmso"],
        help="Columns in metadata to group by for scaling"
    )
    parser.add_argument(
        "--control-label",
        default="DMSO",
        help="Label in 'compound' column to identify controls"
    )

    args = parser.parse_args()

    main(
        metadata_file=args.metadata,
        features_file=args.features,
        output_file=args.output,
        groupby_keys=args.groupby,
        control_label=args.control_label
    )
