import subprocess,os
import pandas as pd
import pysam
from st_cal import cal_base_ratio

def run_base_stat(fasta, outdir,ref, start=None, end=None, vref=False):
    subprocess.run(['samtools', 'faidx', f'{fasta}'], check=True,
                   timeout=3600)  #samtools 建立索引
    pfa = pysam.FastaFile(fasta)
    if vref == False:
        fa_seq = pfa.fetch(reference=ref, start=start, end=end)
        ref_stat = cal_base_ratio(fa_seq)
        AT = ref_stat["AT"]
        GC = ref_stat["GC"]
        pd.DataFrame({
            "base_type": ref_stat.keys(),
            "count_ratio": ref_stat.values()
        }).to_csv(
            os.path.join(outdir, "ref_base_stat.csv"), index=False)
    else:
        references = pfa.references
        references.remove(ref)
        A = T = G = C = other = 0
        for i in references:
            fa_seq = pfa.fetch(i)
            ref_stat_sub = cal_base_ratio(fa_seq)
            A += ref_stat_sub["A"] + A
            T += ref_stat_sub["T"] + T
            G += ref_stat_sub["G"] + G
            C += ref_stat_sub["C"] + C
            other += ref_stat_sub["other"] + other
        total_len = sum(pfa.lengths) - len(pfa.fetch(ref))
        AT = (A + T) / total_len
        GC = (G + C) / total_len
        ref_stat = {
            "A": A,
            "T": T,
            "G": G,
            "C": C,
            "AT": AT,
            "GC": GC,
            "other": other
        }
        pd.DataFrame({
            "base_type": ref_stat.keys(),
            "count_ratio": ref_stat.values()
        }).to_csv(
            os.path.join(outdir, "ref_base_stat.csv"), index=False)
    pfa.close()
    return AT, GC
