import sys 


def usage():
    print('Usage: python script.py [sam file] [outfile_name]')

def calculate_ribosome_density_ontrans(input_file, output_file, min_value=23, max_value=35):
    # save in dict
    density_dict = {}
    total_reads = 0

    # open sam file
    with open(input_file, "r") as reader:
        for line in reader:
            record = line.strip().split("\t")
            # do something
            if record[1] == "0":
                total_reads += 1
                # tags
                refname = record[2]
                align_pos = int(record[3])
                read_length = len(record[9])
            
                # filter read length
                if min_value < read_length < max_value:
                    end5 = align_pos
                    end3 = end5 + read_length - 1
                    # shift +- 11nt
                    centerEnd5 = end5 + 11
                    centerEnd3 = end3 - 11
                    centerLength = centerEnd3 - centerEnd5 + 1
                
                    # ribo density
                    for elem in range(centerEnd5, centerEnd3 + 1):
                        key = f"{refname}:{elem}"
                        density_dict[key] = density_dict.get(key, 1.0 / centerLength) + 1.0 / centerLength
                else:
                    continue
    
    # output file
    with open(output_file, "w") as outfile:
        for key in sorted(density_dict.keys()):
            refname, align_pos = key.split(":")
            raw_density = density_dict[key]
     
            # RPM normalization
            rpm = (raw_density/total_reads) * 1000000
            outfile.write(f"{refname}\t{align_pos}\t{raw_density}\t{rpm}\n")


if __name__ == '__main__':
    try:
        calculate_ribosome_density_ontrans(input_file = sys.argv[1], output_file= sys.argv[2])

    except:
        usage()