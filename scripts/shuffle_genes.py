import pandas as pd
import numpy as np
import sys

# Set seed for reproducibility
np.random.seed(42)

# Read the data from stdin
data = pd.read_table(sys.stdin, index_col=0)

# Get shuffled column indices
shuffled_cols = np.random.permutation(data.columns)

# Shuffle columns
shuffled_data = data[shuffled_cols].copy()

# Reset index to turn sample IDs into a column
shuffled_data = shuffled_data.reset_index()
shuffled_data = shuffled_data.rename(columns={'index': 'Sample'})

# Write the shuffled data to stdout
shuffled_data.to_csv(sys.stdout, sep='\t', index=False)
