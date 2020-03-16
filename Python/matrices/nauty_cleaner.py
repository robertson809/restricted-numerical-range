import sys
name = sys.argv[1]
# name = 'adj6_pt4.txt'

infile = open(name)
outfile = open(name[:-4] + '_cleaned.txt', 'w+')
lineList = infile.readlines()
# outfile.write("%d " % g[i, j])
for row in lineList:
    if row[0] == "G" and len(row) > 1:
        continue
    else:
        outfile.write(row)
outfile.write('\n')
