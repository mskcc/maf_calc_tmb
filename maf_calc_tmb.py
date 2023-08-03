#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to calculate TMB tumor mutational burden value

num_variants / bases_covered * 1,000,000 = TMB in Megabases
"""
import argparse

#using the bed file, function creates the list of covered locations, in this format chr_position
def get_panel_coverage_list(targets_bed_file):
    panel_list=[]
    with open(targets_bed_file) as bed_file:
        for line in bed_file:
            if line.startswith('#'): continue #skip comment lines
            chromosome,spos,epos=line.split()[:3]
            panel_list+= [chromosome+'_'+str(pos) for pos in range(int(spos),int(epos)+1)]           
    return set(panel_list)

def count_on_target_mutations(maf_file,panel_set):
    mutation_set=set()
    f=open(maf_file).read().split('\n')
    for line in f:
        if line.startswith('#'): continue #skip comment lines
        if line.startswith('Hugo_Symbol'):  continue #skip header line
        mutation_set.add('_'.join(line.split()[4:6]))
    return len(set.intersection(mutation_set,panel_set)),len(set.difference(mutation_set,panel_set)),len(mutation_set)

def main():
    """
    Script to calculate TMB tumor mutational burden value
    """
    parser = argparse.ArgumentParser(description='Script to calculate TMB tumor mutational burden value')
    parser.add_argument("--maf_file",         required=True, dest = 'maf_file',         help="Input maf file")
    parser.add_argument("--targets_bed_file", required=True, dest = 'targets_bed_file', help="targets file")
    parser.add_argument("--targets_length",   required=True, dest = 'targets_length',   help="targets_length")
    parser.add_argument("--output_filename",  required=False,dest = 'output_filename',  help="output_filename")
    args = parser.parse_args()

    panel_covered_locations=get_panel_coverage_list(args.targets_bed_file)

    on_target_mutations_count,off_target_mutations_count,total_num_of_mutations=count_on_target_mutations(args.maf_file,panel_covered_locations)    
    
    tmb_val=round(on_target_mutations_count/int(args.targets_length)*1000000,3)
    
    print("Total number of unique mutations:",total_num_of_mutations)
    print("Number of on-target mutations:",on_target_mutations_count)
    print("Number of off-target mutations:",off_target_mutations_count)
    
    print("TMB value:",tmb_val)

    if args.output_filename is not None:
        open(args.output_filename,'w').write(str(tmb_val)+'\n')


if __name__ == '__main__':
    main()