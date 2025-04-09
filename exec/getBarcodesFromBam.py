import pysam
import argparse
import sys
import os

parser = argparse.ArgumentParser(description='Get sam entries for each cell barcode.')
parser.add_argument('--cells',
                    '-c',
                    type=str,
                    default='/home/gdrobertslab/mvc002/analyses/roberts/dev/testSnps/output/two_cancers/two_cancers_clusters_cells_S0149.txt',
                    help='list of cell barcodes with no header')
parser.add_argument('--bam',
                    '-b',
                    type=str,
                    default='/home/gdrobertslab/lab/Counts_2/S0149/possorted_genome_bam.bam',
                    help='bam file output from cellranger')
parser.add_argument('--out_dir',
                    '-o',
                    type=str,
                    default='test',
                    help='output bam file directory')
parser.add_argument('--print_buffer_size',
                    type=int,
                    default=1000,
                    help='Number of sam entries per cell barcode to store before printing them out')
parser.add_argument('--verbose',
                    action='store_true',
                    help='print out extra information')

args = parser.parse_args()

################################################################################
### Global variables
label_dict = {}
umi_dict = {}
print_buffer_dict = {}
bam_outs = {}
################################################################################
### Code

########
### functions

def add_to_label_dict(x):
    cell, label = x.split('\t')
    label_dict[cell] = label
    return(cell)

def populate_write_buffer(cell_barcodes):
    for one_cell_barcode in cell_barcodes:
        print_buffer_dict[one_cell_barcode] = []

def open_bam_outs_from_labels(labels, bam_template):
    for label in labels:
        # check if the output bam file already exists
        bam_file = args.out_dir + '/' + label + '.bam'
        if os.path.exists(bam_file):
            print("Error: bam file already exists", file = sys.stderr)
            sys.exit(1)
        else:
            bam_outs[label] = \
                pysam.AlignmentFile(bam_file, 'wb', template = bam_template)
    return

# For a single bam line, get the CB:Z and UB:Z tags, then if the CB:Z tag is in the
# label_dict as a key, write the line to the appropriate bam file if the UB:Z
# tag has not been seen before (not in the umi_dict dictionary)
def process_line(line):
    # check if CB:Z tag is present
    if line.has_tag('CB') and line.has_tag('UB'):
        # get CB:Z and UB:Z tags
        one_cell_barcode = line.get_tag('CB')
        umi = line.get_tag('UB')
        molecule = \
            one_cell_barcode + \
            umi + \
            str(line.reference_name) + \
            str(line.reference_start)

        # check if CB:Z tag is in label_dict
        if one_cell_barcode in label_dict:
            # check if UB:Z tag is in umi_dict
            if molecule not in umi_dict:
                # add UB:Z tag to umi_dict
                umi_dict[molecule] = 1

                print_buffer_dict[one_cell_barcode].append(line)
                if len(print_buffer_dict[one_cell_barcode]) >= args.print_buffer_size:
                    write_one_cb(one_cell_barcode)

    return

def write_one_cb(one_cell_barcode):
    [
        bam_outs[label_dict[one_cell_barcode]].write(line)
        for line in print_buffer_dict.get(one_cell_barcode)
    ]
    print_buffer_dict[one_cell_barcode] = []

# make main function
def main():
    in_bam = pysam.AlignmentFile(args.bam, 'rb')

    # Read in cell barcodes
    # This assumes there is no header in barcode file
    barcode_file = open(args.cells, 'r')
    cell_barcodes = [add_to_label_dict(x.strip()) for x in barcode_file.readlines()]

    # populate write buffer dict for each cell barcode with an empty list
    populate_write_buffer(cell_barcodes)

    # This gets all unique labels (groups from barcode_file) so we can name one bam file per group
    all_labels = list(set(label_dict.values()))
    # check if cell_barcodes is empty
    if not cell_barcodes:
        print("Error: cell_barcodes is empty", file = sys.stderr)
        sys.exit(1)

    open_bam_outs_from_labels(all_labels, in_bam)

    # Loop through the bam file and use process_line on each line
    # This will write out the data for a cell barcode once there are args.print_buffer_size reads for that cell barcode
    [process_line(line) for line in in_bam]

    # Write out the remaining reads
    [write_one_cb(one_cell_barcode) for one_cell_barcode in cell_barcodes]

if __name__ == '__main__':
    main()