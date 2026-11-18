#!/usr/bin/env python3
"""Regression tests for evidence grading rules, curation fixes, and BMR construct export."""
import csv
import json
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

from build_evidence_six import OUT, grade_record


def read_tsv(path: Path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def main():
    base = {'原证据标签': 'PARTIAL', '支持关系': '实验支持', '功能ID': 'Ni01', '体系ID': 'x', '证据级别': 'C 结构或间接实验'}
    unit_cases = [
        ('纯化蛋白并测定晶体结构', 'E3'),
        ('纯化并测定X射线晶体结构及活性中心配位', 'E3'),
        ('纯化Frh并解析十二聚复合体的冷冻电镜结构，定位F420结合模块', 'E3'),
        ('Fd3与多个Mbh亚基共纯化/质谱检出。', 'E3'),
        ('纯化HoxBC与HoxJ，分析复合体组成、相互作用及H2感知相关信号连接。', 'E3'),
        ('膜复合体活性染色、蛋白鉴定及共营养条件表达', 'E3'),
        ('HoxEFU氧化还原模块研究。', 'E3'),
        ('纯化Hnd与本菌FdxB，在PFOR还原Fd及NADH共同存在时测H2生成；另有ΔhndD及补回。', 'E1'),
        ('野生型/Δech无胞质膜囊泡对照，Fdred驱动H2及ATP生成。', 'E1'),
        ('多组学及E. limosum异源表达；在E. limosum中异源导入pJIR-NADP进行遗传干预，摄气和生长提高。', 'E2'),
        ('比较酸碱条件下hyc、hyf及F0F1相关突变对H2产生的影响。', 'E2'),
    ]
    for text, expected in unit_cases:
        level, _, _ = grade_record({**base, '实际测得': text})
        assert level == expected, (text, level, expected)

    evidence = {r['证据ID']: r for r in read_tsv(OUT / 'evidence.tsv')}
    assert evidence['Ni01-EV-8efaab5cf1']['统一证据等级'] == 'E3'
    assert evidence['Fe08-EV-c7af04ff3c']['统一证据等级'] == 'E2'
    assert evidence['Ni06-EV-d745508085']['统一证据等级'] == 'E1'
    assert evidence['Fe02-EV-2a559302a8']['统一证据等级'] == 'E1'
    assert 'NADPH→CrFNR→CrPetF驱动产氢' not in evidence['Fe02-EV-2a559302a8']['实际测得']
    assert '643/644' in evidence['Ni21-EV-c41776bb3c']['体系序列限制']
    assert 'FRAG_BG04_163_BMR_473_1049' in evidence['Ni21-EV-c41776bb3c']['体系序列限制']
    assert '属水平参考' in evidence['Ni21-EV-71f9a5f6d3']['限制及未测部分']

    summary = {r['功能ID']: r for r in read_tsv(OUT / 'completion_summary.tsv')}
    assert len(summary) == 34
    assert all(int(r['仍需补齐研究数']) == 0 and int(r['六项资料齐全研究数']) >= 6 for r in summary.values())
    assert int(summary['Ni21']['六项资料齐全研究数']) == 6
    assert int(summary['Fe02']['六项资料齐全研究数']) == 7

    frags = read_tsv(OUT / 'construct_fragments.tsv')
    bmr = next(r for r in frags if r['片段ID'] == 'FRAG_BG04_163_BMR_473_1049')
    assert int(bmr['片段长度_aa']) == 577 and int(bmr['片段CDS长度_nt']) == 1731
    assert str(Seq(bmr['片段CDS序列']).translate()) == bmr['片段蛋白序列']
    faa = {r.id: str(r.seq) for r in SeqIO.parse(OUT / 'construct_fragments.faa', 'fasta')}
    fna = {r.id: str(r.seq) for r in SeqIO.parse(OUT / 'construct_fragments_cds.fna', 'fasta')}
    assert faa['FRAG_BG04_163_BMR_473_1049'] == bmr['片段蛋白序列']
    assert fna['FRAG_BG04_163_BMR_473_1049'] == bmr['片段CDS序列']

    print(json.dumps({
        'status': 'passed',
        'unit_regression_cases': len(unit_cases),
        'classes_complete_at_least_6': sum(int(r['六项资料齐全研究数']) >= 6 for r in summary.values()),
        'construct_fragments_exported': len(frags),
    }, ensure_ascii=False))


if __name__ == '__main__':
    main()
