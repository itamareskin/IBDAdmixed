__author__ = 'Itamar'
import os
from itertools import islice

def convert_hapmap_to_ped(hapmap_file, ped_file, map_file):
    with open(hapmap_file) as file_hapmap, open(ped_file, "w") as file_ped, open(map_file, "w") as file_map:
        done = False
        buffer_size = 100
        snp_idx = 0

        line = file_hapmap.readline()
        line = line.split()
        nr_inds = int(len(line)/2-1)
        ids = line[2:]
        rsids = []
        positions = []
        haplos = [[]]*nr_inds

        while True:

            if done:
                break
            # read next buffer_size lines from the file
            lines = list(islice(file_hapmap, buffer_size))

            if len(lines) == 0:
                done = True

            for line in lines:
                #print "line: " + line
                line = line.split()

                rsids.append(line[0])
                positions.append(int(line[1]))
                for ind_idx in range(nr_inds):
                    haplos[ind_idx].append([line[2+2*ind_idx],line[3+2*ind_idx]])

            snp_idx += 1

        nr_snps = snp_idx

        for snp_idx in range(nr_snps):
            file_map.write("1 " + rsids[snp_idx] + " 0 " + str(positions[snp_idx]) + "\n")

        for ind_idx in range(nr_inds):
            file_ped.write(ids[ind_idx] + " " + ids[ind_idx] + " 0 0 1 1")
            for snp_idx in range(nr_snps):
                file_ped.write(" " + haplos[ind_idx][snp_idx][0] + " " + haplos[ind_idx][snp_idx][1])
            file_ped.write("\n")

        x = 1