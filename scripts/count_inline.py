import argparse
import re
import sys

parser = argparse.ArgumentParser(description='Count lines in a file but pass on the data.')
parser.add_argument('-i',
                    '--ignore_regex',
                    type=str,
                    help='Pattern to find lines not to count')
parser.add_argument('-o',
                    '--output',
                    type=str,
                    help='output file',
                    required=True)

args = parser.parse_args()

################################################################################
### Global variables
line_count = 0

################################################################################
### Code

########
# Read in from stdin, print each line to stdout and count the lines
# Only count lines that don't match the ignore_regex
for line in sys.stdin:
    if args.ignore_regex:
        if not re.search(args.ignore_regex, line):
            line_count += 1
    else:
        line_count += 1

    sys.stdout.write(line)

# Write counts to output file
out_file = open(args.output, 'w')
out_file.write(str(line_count))
out_file.close()