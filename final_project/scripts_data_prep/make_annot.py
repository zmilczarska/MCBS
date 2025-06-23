#!/usr/bin/env python
from __future__ import print_function
import pandas as pd
import numpy as np
import argparse
import os
import sys
from pybedtools import BedTool
import gzip

def validate_bed(bed_path):
    """Validate BED file exists and is readable"""
    if not os.path.exists(bed_path):
        raise FileNotFoundError(f"BED file not found: {bed_path}")
    try:
        test = BedTool(bed_path)
        return True
    except Exception as e:
        print(f"Invalid BED file: {e}", file=sys.stderr)
        return False

def gene_set_to_bed(args):
    """Convert gene set to BED format"""
    print('Making gene set bed file')
    GeneSet = pd.read_csv(args.gene_set_file, header=None, names=['GENE'])
    all_genes = pd.read_csv(args.gene_coord_file, delim_whitespace=True)
    
    # Merge and calculate coordinates
    df = pd.merge(GeneSet, all_genes, on='GENE', how='inner')
    df['START'] = np.maximum(1, df['START'] - args.windowsize)
    df['END'] = df['END'] + args.windowsize
    
    # Create BED entries with consistent chromosome formatting
    iter_df = [[f"chr{str(x1)}", x2 - 1, x3] 
              for (x1, x2, x3) in np.array(df[['CHR', 'START', 'END']])]
    
    return BedTool(iter_df).sort().merge()

def make_annot_files(args, bed_for_annot):
    """Create annotation files from BED data"""
    print('Making annot file')
    
    # Read BIM file
    df_bim = pd.read_csv(args.bimfile,
                        delim_whitespace=True,
                        usecols=[0, 1, 2, 3],
                        names=['CHR', 'SNP', 'CM', 'BP'])
    
    # Create BED entries for BIM positions
    iter_bim = [[f"chr{str(x1)}", x2 - 1, x2] 
               for (x1, x2) in np.array(df_bim[['CHR', 'BP']])]
    
    # Find intersections
    bimbed = BedTool(iter_bim)
    annotbed = bimbed.intersect(bed_for_annot)
    bp = [x.start + 1 for x in annotbed]
    
    # Create annotation DataFrame
    df_int = pd.DataFrame({'BP': bp, 'ANNOT': 1})
    df_annot = pd.merge(df_bim, df_int, how='left', on='BP')
    df_annot.fillna(0, inplace=True)
    df_annot = df_annot[['ANNOT']].astype(int)
    
    # Write output
    if args.annot_file.endswith('.gz'):
        with gzip.open(args.annot_file, 'wb') as f:
            df_annot.to_csv(f, sep="\t", index=False)
    else:
        df_annot.to_csv(args.annot_file, sep="\t", index=False)

if __name__ == '__main__':
    # Argument parsing
    parser = argparse.ArgumentParser(description='Create annotation files for LDSC analysis')
    parser.add_argument('--gene-set-file', type=str, 
                       help='File containing gene names (one per line)')
    parser.add_argument('--gene-coord-file', type=str, default='ENSG_coord.txt',
                       help='File with gene coordinates (columns: GENE, CHR, START, END)')
    parser.add_argument('--windowsize', type=int,
                       help='Base pairs to add around transcribed regions')
    parser.add_argument('--bed-file', type=str,
                       help='UCSC BED file with annotation regions')
    parser.add_argument('--nomerge', action='store_true', default=False,
                       help="Don't merge overlapping BED regions")
    parser.add_argument('--bimfile', type=str,
                       help='PLINK BIM file for LD score calculation')
    parser.add_argument('--annot-file', type=str,
                       help='Output annotation file name')
    
    args = parser.parse_args()

    # Input validation
    if args.bed_file and not validate_bed(args.bed_file):
        sys.exit(1)

    try:
        # Process input data
        if args.gene_set_file is not None:
            bed_for_annot = gene_set_to_bed(args)
        else:
            # Try multiple sorting methods
            try:
                # Method 1: Simple sort
                bed_for_annot = BedTool(args.bed_file).sort()
            except Exception as e:
                print(f"Standard sort failed: {e}", file=sys.stderr)
                try:
                    # Method 2: Chromosome-aware sort
                    bed_for_annot = BedTool(args.bed_file).sort(g='hg19.genome')
                except Exception as e:
                    print(f"Chromosome-aware sort failed: {e}", file=sys.stderr)
                    sys.exit(1)
            
            if not args.nomerge:
                bed_for_annot = bed_for_annot.merge()

        # Generate output
        make_annot_files(args, bed_for_annot)
        
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)
