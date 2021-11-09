
import gzip

class Get96rs:
    def __init__(self, invcf, rsfn, outfn, outst):
        self.rs_dic = self.init_rs_dic(rsfn)
        self.add_info(invcf)
        self.write_rs(outfn)
        self.write_stats(outst)

    def write_rs(self, outfn):
        outfh = open(outfn, 'w')
        headers = ['RS']
        headers.append('GT')
        headers.append('AD')
        headers.append('DP')
        headers.append('GQ')
        headers.append('PL')
        outfh.write('{0}\n'.format('\t'.join(headers)))
        for rs, info_dic in self.rs_dic.items():
            items = [rs]
            items.append(info_dic['GT'])
            items.append(info_dic['AD'])
            items.append(info_dic['DP'])
            items.append(info_dic['GQ'])
            items.append(info_dic['PL'])
            outfh.write('{0}\n'.format('\t'.join(items)))
        outfh.close()

    def write_stats(self, outst):
        call = 0
        miss = 0
        for rs, info_dic in self.rs_dic.items():
            if info_dic['GT'] in ['NoCall']:
                miss += 1
            else:
                call += 1
        outfh = open(outst, 'w')
        outfh.write(f'Call : {call}\n')
        outfh.write(f'Miss : {miss}\n')
        outfh.close()

    def add_info(self, invcf):
        for line in gzip.open(invcf, 'rt'):
            if line.startswith('##'):
                continue
            items = line.rstrip('\n').split('\t')
            if items[0] in ['#CHROM']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            rs = items[idx_dic['ID']]
            if rs in self.rs_dic:
                ref = items[idx_dic['REF']]
                alt = items[idx_dic['ALT']]
                info_dic = self.get_info(items[idx_dic['FORMAT']], items[9])
                gt = self.get_gt(ref, alt, info_dic)
                for fmt, value in info_dic.items():
                    if fmt in ['GT']:
                        self.rs_dic[rs][fmt] = gt
                    else:
                        self.rs_dic[rs][fmt] = value
            else:
                continue
        return 1

    def get_gt(self, ref, alt, info_dic):
        gt_s = [ref]
        gt_s.extend(alt.split(','))

        gt_idx = info_dic['GT']
        gt_idx = gt_idx.replace('|','/')
        gt_idx_s = gt_idx.split('/')

        gt_str = list()
        for gt_idx in [x for x in gt_idx_s]:
            if gt_idx in ["."]:
                gt_str.append("N")
            else:
                gt_str.append(gt_s[int(gt_idx)])
        return '/'.join(gt_str)

    def get_info(self, _fmt, _info):
        infos = _info.split(':')
        info_dic = dict()
        for idx, fmt in enumerate(_fmt.split(':')):
            info_dic.setdefault(fmt, infos[idx])
        return info_dic


    def init_rs_dic(self, rsfn):
        rs_dic = dict()
        for line in open(rsfn):
            items = line.rstrip('\n').split('\t')
            rs_dic.setdefault(items[0], {}).setdefault('GT', 'NoCall')
            rs_dic.setdefault(items[0], {}).setdefault('AD', 'NoCall')
            rs_dic.setdefault(items[0], {}).setdefault('DP', 'NoCall')
            rs_dic.setdefault(items[0], {}).setdefault('GQ', 'NoCall')
            rs_dic.setdefault(items[0], {}).setdefault('PL', 'NoCall')
        return rs_dic


def main(args):
    get96rs = Get96rs(snakemake.input["vcf"],
                      snakemake.params["rsfn"],
                      snakemake.output["rs"],
                      snakemake.output["stats"])

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--invcf')
    parser.add_argument('--rsfn')
    parser.add_argument('--outfn')
    args = parser.parse_args()
    main(args)
