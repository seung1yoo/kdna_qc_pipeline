


import os
import yaml

class MergeResults:
    def __init__(self, args):
        self.rs_s = list()
        self.r_dic = self.init_r_dic(args.config)
        for sample, info_dic in self.r_dic.items():
            self.add_bases_mapped(args.result_path, sample)
            self.add_96_stats(args.result_path, sample)
            self.add_96_genotype(args.result_path, sample)
        self.write_result(args.outprefix)

    def init_r_dic(self, config):
        r_dic = dict()
        c_dic = yaml.load(open(config))
        for sample, lane_dic in c_dic['samples'].items():
            r_dic.setdefault(sample, {}).setdefault('bases mapped', 0)
            r_dic.setdefault(sample, {}).setdefault('call', 0)
            r_dic.setdefault(sample, {}).setdefault('miss', 0)
            r_dic.setdefault(sample, {}).setdefault('rs', {})
        return r_dic

    def add_bases_mapped(self, r_path, sample):
        for line in open(os.path.join(r_path, sample, f"{sample}.align.stats")):
            if line.startswith('bases mapped:'):
                #bases mapped:   105565005698    # ignores clipping
                value = int(line.split()[2])
                self.r_dic[sample]['bases mapped'] = value
            else:
                continue
        return 1

    def add_96_stats(self, r_path, sample):
        for line in open(os.path.join(r_path, sample, f"{sample}.96.stats")):
            items = line.rstrip('\n').split(':')
            items = [item.strip() for item in items]
            if items[0] in ['Call']:
                self.r_dic[sample]['call'] = int(items[1])
            elif items[0] in ['Miss']:
                self.r_dic[sample]['miss'] = int(items[1])
            else:
                continue
        return 1

    def add_96_genotype(self, r_path, sample):
        for line in open(os.path.join(r_path, sample, f"{sample}.96.rs")):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['RS']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            rs = items[idx_dic['RS']]
            if rs not in self.rs_s:
                self.rs_s.append(rs)
            gt = items[idx_dic['GT']]
            ad = items[idx_dic['AD']]
            dp = items[idx_dic['DP']]
            gq = items[idx_dic['GQ']]
            pl = items[idx_dic['PL']]
            self.r_dic[sample]['rs'].setdefault(rs, {}).setdefault('gt', gt)
            self.r_dic[sample]['rs'].setdefault(rs, {}).setdefault('ad', ad)
            self.r_dic[sample]['rs'].setdefault(rs, {}).setdefault('dp', dp)
            self.r_dic[sample]['rs'].setdefault(rs, {}).setdefault('gq', gq)
            self.r_dic[sample]['rs'].setdefault(rs, {}).setdefault('pl', pl)
        return 1

    def write_result(self, outprefix):
        outfh = open(f"{outprefix}.qc.tsv", "w")
        outfh.write('SampleName\tBasesMapped\tVarCall\tVarMiss\n')
        for sample, info_dic in self.r_dic.items():
            items = [sample]
            items.append(info_dic['bases mapped'])
            items.append(info_dic['call'])
            items.append(info_dic['miss'])
            outfh.write("{0}\n".format('\t'.join([str(x) for x in items])))
        outfh.close()

        outfh = open(f"{outprefix}.rs.tsv", "w")
        items = ['SampleName']
        for rs in self.rs_s:
            for info in ['gt','ad','dp','gq','pl']:
                items.append(f"{rs}.{info}")
        outfh.write("{0}\n".format('\t'.join(items)))
        for sample, info_dic in self.r_dic.items():
            items = [sample]
            for rs in self.rs_s:
                for info in ['gt','ad','dp','gq','pl']:
                    items.append(info_dic['rs'][rs][info])
            outfh.write("{0}\n".format('\t'.join([str(x) for x in items])))
        outfh.close()

        outfh = open(f"{outprefix}.rs.gt.tsv", "w")
        items = ['SampleName']
        items.extend(self.rs_s)
        outfh.write("{0}\n".format('\t'.join(items)))
        for sample, info_dic in self.r_dic.items():
            items = [sample]
            for rs in self.rs_s:
                items.append(info_dic['rs'][rs]['gt'])
            outfh.write("{0}\n".format('\t'.join([str(x) for x in items])))
        outfh.close()





def main(args):

    result_obj = MergeResults(args)




if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--result-path', default='/data08/project/dev_test_overload/analysis/result')
    parser.add_argument('-c', '--config', default='/data08/project/dev_test_overload/dev.config.yaml')
    parser.add_argument('-o', '--outprefix', default='test')
    args = parser.parse_args()
    main(args)
