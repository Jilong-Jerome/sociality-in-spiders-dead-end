import pandas as pd
import numpy as np
import sys

# Custom interpolation function
def custom_interpolate(df, column):
    filled_values = []
    prev_row = None
    next_row = None
    
    for _, row in df.iterrows():
        if pd.isna(row[column]):
            if prev_row is not None and next_row is None:
                next_row_candidate = df.loc[(df[column].notna()) & (df['years_before_present'] > row['years_before_present'])]
                if not next_row_candidate.empty:
                    next_row = next_row_candidate.iloc[0]
            if prev_row is not None and next_row is not None:
                time_diff = next_row['years_before_present'] - prev_row['years_before_present']
                value_diff = next_row[column] - prev_row[column]
                time_fraction = (row['years_before_present'] - prev_row['years_before_present']) / time_diff
                filled_values.append(prev_row[column] + value_diff * time_fraction)
                if row['years_before_present'] > next_row['years_before_present']:
                    prev_row = next_row
                    next_row = None
            else:
                filled_values.append(0)
        else:
            filled_values.append(row[column])
            prev_row = row
            next_row = None
            
    return filled_values

# Read the input TSV files
autosomal_file = sys.argv[1]
x_chromosome_file = sys.argv[2]

autosomal_df = pd.read_csv(autosomal_file, sep='\s+', names=['years_before_present', 'Ne_Auto'], usecols=[0, 1])
x_chromosome_df = pd.read_csv(x_chromosome_file, sep='\s+', names=['years_before_present', 'Ne_X'], usecols=[0, 1])

# Merge the dataframes on 'years_before_present' and sort them
merged_df = autosomal_df.merge(x_chromosome_df, on='years_before_present', how='outer').sort_values('years_before_present').reset_index(drop=True)

# Fill the missing values using the custom interpolation function
merged_df['Ne_Auto'] = custom_interpolate(merged_df, 'Ne_Auto')
merged_df['Ne_X'] = custom_interpolate(merged_df, 'Ne_X')

# Replace any remaining NaN values with 0
merged_df.fillna(0, inplace=True)

# Calculate the ratio between Ne_X and Ne_Auto
merged_df['Ne_X/Ne_Auto'] = merged_df['Ne_X'] / merged_df['Ne_Auto']

# Add the 'species' column and fill it with the given user input
species_name = sys.argv[3]
merged_df['species'] = species_name

# Save the merged dataframe as a new TSV file
output_file = sys.argv[3]+"_correct.tsv"
merged_df.to_csv(output_file, sep='\t', index=False)

print("Merged file saved to", output_file)

