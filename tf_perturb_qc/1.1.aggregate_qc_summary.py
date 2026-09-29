#!/usr/bin/env python
# %%
import glob
import os

import pandas as pd

# %%
core_path = "/oak/stanford/groups/engreitz/Users/opushkar/tf_perturb"
output_path = os.path.join(core_path, "2.qc_gex/1.qc_data")

days = ["d0", "d1", "d2", "d3"]
reps = [1, 2]

# %%
expected_tags = [f"{day}_rep{rep}" for day in days for rep in reps]
per_tag_files = {
    tag: os.path.join(output_path, f"qc_filtering_summary_{tag}.tsv") for tag in expected_tags
}
missing = [tag for tag, path in per_tag_files.items() if not os.path.exists(path)]
assert not missing, f"Missing per-tag QC summaries, array task(s) did not complete: {missing}"

summary_df = pd.concat([pd.read_csv(path, sep="\t") for path in per_tag_files.values()], ignore_index=True)
summary_df = summary_df.sort_values(["day", "rep"]).reset_index(drop=True)
summary_df.to_csv(os.path.join(output_path, "qc_filtering_summary.tsv"), sep="\t", index=False)

print("QC filtering summary:")
print(summary_df.to_string(index=False))
# %%
