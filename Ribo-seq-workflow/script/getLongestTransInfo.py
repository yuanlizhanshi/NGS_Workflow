import gzip
import sys
# define function
def usage():
    print('Usage: python script.py [gtf file]')

def getLongestTransInfo(gtf_file,longest_file):
    #####################################################################################
    # 1 extract transcript features information
    #####################################################################################
    # parse gtf file
    if gtf_file.endswith(".gz"):
        gtf =  gzip.open(gtf_file,"rt")
    else:
        gtf =  open(gtf_file,"r")
        
    trans_5utr = {}
    trans_cds = {}
    trans_3utr = {}
    trans_exon = {}
    gene_with_trans = {}
        
    # loops
    for line in gtf:
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        # attributes for gtf
        attributes = fields[8]
        # get feature type
        feature = fields[2]
        
        # transcript id with cds length
        # _,gene_id,transcript_id = extract_record(attributes)
        try:
            transcript_id = attributes.split('transcript_id "')[1].split('"')[0]
        except IndexError:
            transcript_id = "noTransId"
        f_len = abs(int(fields[4]) - int(fields[3])) + 1
        
        # get type of line
        if feature == "CDS":
            # gene_id with transcript
            try:
                gene_id = attributes.split('gene_id "')[1].split('"')[0]
            except IndexError:
                gene_id = "noGeneId"
                
            if gene_id not in gene_with_trans:
                gene_with_trans[gene_id] = [transcript_id]
            else:
                if transcript_id not in gene_with_trans[gene_id]:
                    gene_with_trans[gene_id].append(transcript_id)
                
            trans_cds.setdefault(transcript_id,0)
            trans_cds[transcript_id] += f_len
        # save 5utr/3utr length
        elif feature in ("5UTR","five_prime_utr"):
            trans_5utr.setdefault(transcript_id,0)
            trans_5utr[transcript_id] += f_len
        elif feature in ("3UTR","three_prime_utr"):
            trans_3utr.setdefault(transcript_id,0)
            trans_3utr[transcript_id] += f_len
        elif feature == "exon":
            trans_exon.setdefault(transcript_id,0)
            trans_exon[transcript_id] += f_len
        else:
            continue
        
    #####################################################################################
    # 2 extract longest(acorroding to CDS and exon) transcript
    #####################################################################################
    lonest_transcript = []
    for g,trans in gene_with_trans.items():
        # get mached length information
        cds_len = []
        for ti in trans:
            cds_len.append(trans_cds.get(ti,0))
        exon_len = []
        for ti in trans:
            exon_len.append(trans_exon.get(ti,0))
        
        # get max cds length index
        if len(cds_len) > 0 and len(exon_len) > 0:
            # check mutiple equeal cds length
            if len(set(cds_len)) == len(cds_len):
                max_index = cds_len.index(max(cds_len))
                lonest_transcript.append(trans[max_index])
            else:
                # select longest cds length and longest exon length
                max_cds = max(cds_len)
                max_cds_exon = []
                for index,val in enumerate(cds_len):
                    if val == max_cds:
                        max_cds_exon.append(exon_len[index])
                max_exon_index = exon_len.index(max(max_cds_exon))
                lonest_transcript.append(trans[max_exon_index])
        else:
            continue
        
    #####################################################################################
    # 3 extract longest transcript CDS and exon positions
    #####################################################################################
    # parse gtf file
    if gtf_file.endswith(".gz"):
        gtf =  gzip.open(gtf_file,"rt")
    else:
        gtf =  open(gtf_file,"r")
        
    lonest_transcript_set = set(lonest_transcript)
    trans_exon_pos = {}
    trans_cds_pos = {}
    trans_chrom_strand = {}
    
    for line in gtf:
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        # attributes for gtf
        attributes = fields[8]
        # get feature type
        feature = fields[2]
        # transcript id with cds length
        # _, _, transcript_id = extract_record(attributes)
        try:
            transcript_id = attributes.split('transcript_id "')[1].split('"')[0]
        except IndexError:
            continue
        
        # extract transcript exon position and CDS position
        if transcript_id in lonest_transcript_set and feature in ("exon", "CDS"):
            val = f"{fields[3]}:{fields[4]}"
            if feature == "exon":
                trans_exon_pos.setdefault(transcript_id, []).append(val)
            elif feature == "CDS":
                trans_cds_pos.setdefault(transcript_id, []).append(val)
            else:
                continue
        
        # save transcript id chromosome and strand
        if transcript_id in lonest_transcript_set:
            if transcript_id not in trans_chrom_strand:
                try:
                    gene_name = attributes.split('gene_name "')[1].split('"')[0]
                except IndexError:
                    gene_name = "noGeneName"
                try:
                    gene_id = attributes.split('gene_id "')[1].split('"')[0]
                except IndexError:
                    gene_id = "noGeneId"
                chrom = fields[0];strand = fields[6]
                trans_chrom_strand[transcript_id] = [gene_name,gene_id,transcript_id,chrom,strand]
            else:
                continue
        else:
            continue
        
    #####################################################################################
    # 4 output file
    #####################################################################################
    # output save in another file
    output_file = open(longest_file,'w')

    # separate sequences
    for tid,val in trans_chrom_strand.items():
        key_id = "|".join(val)
        val_id = "\t".join(val)
        
        # revese - strand position
        strand = val[4]
        if strand == "-":
            exon_pos = ",".join(trans_exon_pos.get(tid,[])[::-1])
            cds_pos = ",".join(trans_cds_pos.get(tid,[])[::-1])
        else:
            exon_pos = ",".join(trans_exon_pos.get(tid,[]))
            cds_pos = ",".join(trans_cds_pos.get(tid,[]))
        
        line = "\t".join([key_id,val_id,
                        cds_pos,exon_pos,
                        str(trans_5utr.get(tid,0)),str(trans_cds.get(tid,0) + 3),str(trans_3utr.get(tid,0))])
        output_file.write(line + '\n')
        
    # file close
    output_file.close()
    gtf.close()
    
if __name__ == '__main__':
    try:
        getLongestTransInfo(sys.argv[1], 'longest_info.txt')
    except:
        usage()