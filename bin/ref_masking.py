
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def target_pos_finder(bed):
    tpos_dic = dict()
    for line in open(bed):
        if line.startswith('#'):
            continue
        items = line.rstrip('\n').split('\t')
        chrom = items[0]
        start = int(items[1])
        end = int(items[2])
        anno = items[3]
        for pos in range(start, end+1, 1):
            tpos_dic.setdefault(chrom, {}).setdefault(pos, anno)
    return tpos_dic



def main(args):
    tpos_dic = target_pos_finder(args.target)
    for chrom, pos_dic in tpos_dic.items():
        print(chrom, len(pos_dic))

    outfh = open(args.ofa, 'w')
    for record in SeqIO.parse(open(args.ifa), 'fasta'):
        chrom = record.id
        print(chrom)
        if chrom not in tpos_dic:
            continue
        seq = list()
        for idx, nucl in enumerate(str(record.seq)):
            pos = idx+1
            if pos in tpos_dic[chrom]:
                seq.append(nucl)
            else:
                seq.append('N')
        m_record = SeqRecord(Seq(''.join(seq)), id=chrom, description='N_MASKING_SIYOO')
        SeqIO.write(m_record, outfh, 'fasta')
    outfh.close()









if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--ifa')
    parser.add_argument('--ofa')
    parser.add_argument('--target')
    args = parser.parse_args()
    main(args)
