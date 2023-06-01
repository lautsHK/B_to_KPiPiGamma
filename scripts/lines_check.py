import sys

if len(sys.argv) != 2:
    print("Usage: python lines_check.py <filename>")
    sys.exit(1)

filename = sys.argv[1]

with open(filename, 'r') as f:
    for i, line in enumerate(f):
        values = line.strip().split()
        if len(values) != 2:
            print "Line %d contains %d values instead of 2." % (i+1, len(values))
        else:
            for value in values:
                decimal_count = value.count('.')
                if decimal_count != 1:
                    print "Line %d contains value %s with %d decimal points instead of 1." % (i+1, value, decimal_count)
