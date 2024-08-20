#!/usr/bin/python3.8
import mrcfile
import numpy as np
import argparse

def compute_cross_correlation(x1, x2, mask=None):
    if mask is not None:
        mask = mask >= 1  # Keep only values where the mask is >= 1
        x1 = x1[mask]
        x2 = x2[mask]
    
    numerator = np.sum(x1 * x2)
    denominator = np.sqrt(np.sum(x1**2) * np.sum(x2**2))
    return numerator / denominator

def compute_cross_correlation_files(file1, file2, mask_file=None):
    # Load the first MRC file
    with mrcfile.open(file1, mode='r') as mrc1:
        volume1 = mrc1.data

    # Load the second MRC file
    with mrcfile.open(file2, mode='r') as mrc2:
        volume2 = mrc2.data

    # Optionally load the mask file
    if mask_file:
        with mrcfile.open(mask_file, mode='r') as mrc_mask:
            mask = mrc_mask.data
    else:
        mask = None

    # Convert the volumes to numpy arrays and subtract the mean
    x1 = np.array(volume1, dtype=np.float32)
    mean_x1 = np.mean(x1[mask >= 1]) if mask is not None else np.mean(x1)
    x1_centered = x1 - mean_x1

    x2 = np.array(volume2, dtype=np.float32)
    mean_x2 = np.mean(x2[mask >= 1]) if mask is not None else np.mean(x2)
    x2_centered = x2 - mean_x2

    # Compute and return the cross-correlation, considering the mask
    return compute_cross_correlation(x1_centered, x2_centered, mask)

def main():
    # Set up argument parser with detailed description and usage
    parser = argparse.ArgumentParser(
        description='Compute cross-correlation between MRC files with optional masking.',
        usage='%(prog)s [-h] -l LIST -r REFERENCE [-m MASK]',
        epilog='Example usage: python3.8 script.py -l file_list.txt -r reference.mrc -m mask.mrc'
    )
    parser.add_argument('-l', '--list', required=True, help='Path to the text file with the list of MRC files')
    parser.add_argument('-r', '--reference', required=True, help='Path to the reference MRC file')
    parser.add_argument('-m', '--mask', required=False, help='Path to the mask MRC file (optional)')

    args = parser.parse_args()

    # Define the reference file
    reference_file = args.reference

    max_correlation = -1  # Initialize max correlation to a very low value
    max_correlation_file = None
    perfect_matches = []

    # Read the list of files from the text file
    with open(args.list, 'r') as file_list:
        file_paths = [line.strip() for line in file_list]

    # Iterate over the file paths and compare those ending with 'rec.mrc' against the reference
    for file_path in file_paths:
        if file_path.endswith('.mrc'):
            # Compute the cross-correlation with the reference volume, optionally using a mask
            cross_correlation = compute_cross_correlation_files(file_path, reference_file, args.mask)

            # Check if this is the highest correlation we've seen
            if cross_correlation > max_correlation:
                max_correlation = cross_correlation
                max_correlation_file = file_path

            # Check for perfect match (correlation of 1)
            if cross_correlation == 1:
                perfect_matches.append(file_path)

            # Print the result for each file
            print(f"Cross-correlation between {file_path} and {reference_file}: {cross_correlation}")

    # Output the file with the highest correlation
    print(f"\nFile with the highest correlation: {max_correlation_file} (Correlation: {max_correlation})")

    # Output all files with a perfect correlation of 1
    if perfect_matches:
        print("Files with a perfect correlation of 1:")
        for match in perfect_matches:
            print(f" - {match}")
    else:
        print("No files with a perfect correlation of 1.")

if __name__ == '__main__':
    main()
