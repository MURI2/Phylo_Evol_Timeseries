import numpy
import sys
import bz2

input_filename = sys.argv[1]
depth_filename = sys.argv[2]
snp_filename = sys.argv[3]

input_file = bz2.BZ2File(input_filename,"r")
snp_file = bz2.BZ2File(snp_filename,"w")
depth_file = bz2.BZ2File(depth_filename,"w")

avg_depths = None
times = None
alts = None

times_to_change = ['D100', 'D200', 'D300']

depth_records = []

for line in input_file:
    items = line.split(",")
    position = long(items[1])
    allele = items[2].strip()

    if allele[1:3]!='->':
        continue # not a snp!

    for t in times_to_change:
        line = line.replace(t, t[1:])
    snp_file.write(line)

    # calculate depths and add them
    #print [float(subitem[1:]) for subitem in items[3].split()]
    times = numpy.array([float(subitem[1:]) for subitem in items[3].split()])
    depths = [float(subitem) for subitem in items[5].split()]

    depth_records.append(depths)

depths = numpy.array(depth_records)
# Could do median or mean
#avg_depths = depths.mean(axis=0)
avg_depths = numpy.median(depths, axis=0)

alts = numpy.array([0 for t in times])

# 'REL606'
depth_line = ", ".join(["REL606", "0", "Depth", " ".join([str(t) for t in times]), " ".join([str(alt) for alt in alts]), " ".join([str(avg_depth) for avg_depth in avg_depths])])
depth_file.write(depth_line)
depth_file.write("\n")

input_file.close()
snp_file.close()
depth_file.close()
