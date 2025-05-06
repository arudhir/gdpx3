# GDPX3

Right now mostly just me playing around. The released median-collapsed features and metadata can be downloaded from [here](https://datapoints.ginkgo.bio/functional-genomics/gdpx3)

`data/gdpx3_medians.h5ad` is post-processed according to their protocol:

    - Standardized with respect to container and control DMSO
    - Feature selection with PyCytominer

```bash
$ uv run src/process_gdpx3.py \
    --metadata "data/ginkgo-datapoints_gdpx3_data-package/Ginkgo Datapoints_GDPx3_metadata.csv" \
    --features "data/median_features.csv" \
    --output "data/gdpx3_medians.h5ad"

📥 Loading metadata and feature tables...
🔎 Metadata shape: (2528, 15)
🔎 Features shape: (2479, 2769)
🔗 Aligned data on 2,479 sample_ids
📦 Constructing AnnData object...
🧪 AnnData shape: (2479, 2769)
⚖️ Scaling features per group using RobustScaler...
✅ Scaling complete. Scaled data stored in `adata.layers['scaled']`.
🧹 Running feature selection via pycytominer...
✅ Feature selection complete. Final shape: (2479, 766)
🧼 Cleaning obs metadata to ensure h5ad compatibility...
✅ Metadata cleaned.
💾 Writing processed AnnData to data/gdpx3_medians.h5ad...
🎉 Done.
```
