"""
非天然氨基酸（NNAA）文献检索关键词库

设计原则：
1. 按结构母体分组，构建 OR 检索块，再与各轨道工艺词组合
2. 每个化合物收录主名 + 所有常见别称/同义词（IUPAC 名、俗名、缩写、立体前缀变体）
3. 包含 GCE（遗传密码扩展）专属轨道，覆盖 aaRS/tRNA 相关高频 ncAA 化合物
4. 化合物来源：PDF《Derivative of Natural AA》+ iNCLusive 数据库高频化合物
"""

from __future__ import annotations

from itertools import product
from typing import Dict, List, Optional, Sequence

# ===========================================================================
# 第一部分：PDF 化合物名称 → 可检索英文全称 + 全部别称
# （按母体氨基酸分组，非天然/半天然化合物）
# ===========================================================================

# ---------------------------------------------------------------------------
# Lys 衍生类
# ---------------------------------------------------------------------------
LYS_COMPOUND_NAMES: Sequence[str] = (
    # D-Lys
    "D-lysine",
    # Lys(Ac) — N-epsilon-acetyllysine
    "N-epsilon-acetyllysine", "Nepsilon-acetyllysine", "N6-acetyllysine",
    # Lys(Me2) — dimethyllysine
    "N,N-dimethyllysine", "dimethyllysine", "N6,N6-dimethyllysine",
    # Lys(TFA)
    "trifluoroacetyllysine", "N-trifluoroacetyllysine",
    # Dap — 2,3-diaminopropionic acid
    "2,3-diaminopropionic acid", "diaminopropionic acid", "L-2,3-diaminopropionic acid",
    "alpha,beta-diaminopropionic acid", "Dap",
    # Dab — 2,4-diaminobutyric acid
    "2,4-diaminobutyric acid", "diaminobutyric acid", "L-2,4-diaminobutyric acid",
    "alpha,gamma-diaminobutyric acid", "Dab",
    # Orn
    "ornithine", "L-ornithine", "2,5-diaminopentanoic acid",
    # others
    "azetidylglycine", "N-azetidylglycine",
    "morpholinylalanine", "morpholinyl alanine",
)

# ---------------------------------------------------------------------------
# Arg 衍生类
# ---------------------------------------------------------------------------
ARG_COMPOUND_NAMES: Sequence[str] = (
    # hArg — homoarginine
    "homoarginine", "L-homoarginine", "N6-amidinolysin",
    "N6-(aminoiminomethyl)lysine",
    # Cit — citrulline
    "citrulline", "L-citrulline", "2-amino-5-(carbamoylamino)pentanoic acid",
    # h-Cit — homocitrulline
    "homocitrulline", "L-homocitrulline",
    # D-Arg
    "D-arginine",
    # D-norArg, D-HomoArg
    "norarginine", "D-norarginine", "D-homoarginine",
    # Arg(Di-Me) — ADMA
    "asymmetric dimethylarginine", "ADMA",
    "NG,NG-dimethylarginine", "Nomega,Nomega-dimethylarginine",
    # Arg(Me) — MMA
    "monomethylarginine", "methylarginine", "NG-monomethylarginine", "MMA",
    "NG-methylarginine",
    # Arg(NO2) — nitroarginine
    "nitroarginine", "NG-nitro-L-arginine", "L-NAME precursor",
    # Arg(Me2) — SDMA
    "symmetric dimethylarginine", "SDMA",
    "NG,N'G-dimethylarginine", "Nomega,N'omega-dimethylarginine",
    # Arg(c-Pr)
    "cyclopropylarginine",
    # Agp, Agb
    "guanidinopropionic acid", "beta-guanidinopropionic acid",
    "guanidinobutyric acid", "gamma-guanidinobutyric acid",
    # Alb — agmatine-like / 2-amino-4-guanidinobutanoic acid
    "agmatine", "2-amino-4-guanidinobutanoic acid",
    # Phe(4-guan), Phg(4-guan)
    "4-guanidinophenylalanine", "para-guanidinophenylalanine",
    "4-guanidinophenylglycine",
)

# ---------------------------------------------------------------------------
# 脂肪链类（Ala/Val/Leu/Met 骨架衍生）
# ---------------------------------------------------------------------------
ALIPHATIC_COMPOUND_NAMES: Sequence[str] = (
    # D-AA aliphatic series（不对称合成/酶法拆分 D-型氨基酸）
    "D-phenylalanine", "D-leucine", "D-valine", "D-alanine",
    "D-methionine", "D-serine", "D-tyrosine", "D-threonine",
    "D-tryptophan", "D-proline",
    "D-norvaline", "D-norleucine",
    # Nva — norvaline
    "norvaline", "L-norvaline", "2-aminopentanoic acid",
    "alpha-aminovaleric acid", "alpha-aminopentanoic acid",
    # Nle — norleucine
    "norleucine", "L-norleucine", "2-aminohexanoic acid",
    "alpha-aminocaproic acid",
    # Tle — tert-leucine
    "tert-leucine", "L-tert-leucine", "tert-butylglycine",
    "beta,beta-dimethylnorvaline",
    "(2S)-2-amino-3,3-dimethylbutanoic acid",
    "3,3-dimethyl-L-alanine",
    # Aib — alpha-aminoisobutyric acid
    "Aib", "alpha-aminoisobutyric acid", "2-aminoisobutyric acid",
    "alpha-methylalanine", "alpha-methyl alanine",
    "2-amino-2-methylpropanoic acid", "2-amino-2-methylpropionic acid",
    "2-methylalanine",
    # Abu — 2-aminobutyric acid
    "2-aminobutyric acid", "alpha-aminobutyric acid",
    "L-2-aminobutyric acid", "Abu",
    # hLeu — homoleucine
    "homoleucine", "L-homoleucine", "2-amino-4-methylpentanoic acid",
    # NptGly — neopentylglycine
    "neopentylglycine", "N-neopentylglycine",
    # Ahp — 2-aminoheptanoic acid
    "2-aminoheptanoic acid", "alpha-aminoheptanoic acid",
    # Aoc — 2-aminooctanoic acid
    "2-aminooctanoic acid", "aminooctanoic acid",
    # Met(O2) — methionine sulfone
    "methionine sulfone", "L-methionine sulfone", "methionine-S,S-dioxide",
    # Ala(CN) — beta-cyanoalanine
    "beta-cyanoalanine", "3-cyanoalanine", "2-amino-3-cyanopropanoic acid",
    # Algly — allylglycine
    "allylglycine", "L-allylglycine", "2-amino-4-pentenoic acid",
    "alpha-allylglycine",
    # PreGly / Pra — propargylglycine
    "propargylglycine", "L-propargylglycine",
    "2-amino-4-pentynoic acid", "alpha-propargylglycine",
    # MethaGly — methallylglycine
    "methallylglycine", "2-amino-4-methylpent-4-enoic acid",
    # Nva(3-Et) — 3-ethylnorvaline
    "3-ethylnorvaline",
    # Nle(6-OH) — 6-hydroxynorleucine
    "6-hydroxynorleucine",
    # Val(3-OH) / beta-hydroxyvaline
    "3-hydroxyvaline", "beta-hydroxyvaline", "L-beta-hydroxyvaline",
    # homoserine
    "homoserine", "L-homoserine", "2-amino-4-hydroxybutyric acid",
    "DL-homoserine",
    # hSer derivatives
    "O-methylhomoserine", "O-ethylhomoserine",
    # Abu(CF3)
    "trifluoromethylaminobutyric acid", "4,4,4-trifluoro-2-aminobutanoic acid",
)

# ---------------------------------------------------------------------------
# 脂环类（cyclopropyl / cyclobutyl / cyclopentyl / cyclohexyl 骨架）
# ---------------------------------------------------------------------------
CYCLIC_COMPOUND_NAMES: Sequence[str] = (
    # Cha — cyclohexylalanine
    "cyclohexylalanine", "beta-cyclohexylalanine", "L-cyclohexylalanine",
    "3-cyclohexylalanine", "Cha",
    "2-amino-3-cyclohexylpropionic acid",
    # h-Cha — homocyclohexylalanine
    "homocyclohexylalanine", "4-cyclohexyl-2-aminobutanoic acid",
    # cPrA — cyclopropylalanine
    "cyclopropylalanine", "beta-cyclopropylalanine",
    "2-amino-3-cyclopropylpropionic acid",
    # Ala(cBu) — cyclobutylalanine
    "cyclobutylalanine", "2-amino-3-cyclobutylpropanoic acid",
    # cPenA — cyclopentylalanine
    "cyclopentylalanine", "2-amino-3-cyclopentylpropanoic acid",
    # Phg — phenylglycine
    "phenylglycine", "L-phenylglycine", "D-phenylglycine",
    "alpha-phenylglycine", "2-phenylglycine",
    "aminophenylacetic acid", "2-amino-2-phenylacetic acid",
    # Chg — cyclohexylglycine
    "cyclohexylglycine", "alpha-cyclohexylglycine", "Chg",
    "2-amino-2-cyclohexylacetic acid",
    # Gly(cPent) — cyclopentylglycine
    "cyclopentylglycine", "2-amino-2-cyclopentylacetic acid",
    # Gly(cBu) — cyclobutylglycine
    "cyclobutylglycine", "2-amino-2-cyclobutylacetic acid",
    # Igl — indanylglycine
    "indanylglycine", "alpha-indanylglycine",
    "2-(indan-1-yl)glycine", "2-(indan-2-yl)glycine",
    # Aib (also in aliphatic)
    "alpha-aminoisobutyric acid", "alpha-methylalanine", "Aib",
    # ACPrC — 1-aminocyclopropane-1-carboxylic acid
    "1-aminocyclopropane-1-carboxylic acid",
    "1-amino-1-cyclopropanecarboxylic acid",
    "ACC", "cyclopropane amino acid",
    # cHex — 1-aminocyclohexane-1-carboxylic acid
    "1-aminocyclohexane-1-carboxylic acid",
    "aminocyclohexane carboxylic acid", "alpha-aminocyclohexane carboxylic acid",
    # Gly(BCP) — BCP glycine
    "bicyclo[1.1.1]pentane glycine", "BCP glycine", "bicyclopentylglycine",
    "1-(aminomethyl)bicyclo[1.1.1]pentane-1-carboxylic acid",
    # 2F-Cha — difluorocyclohexylalanine
    "difluorocyclohexylalanine", "4,4-difluorocyclohexylalanine",
)

# ---------------------------------------------------------------------------
# 芳香类（Phe / Tyr / His 骨架衍生）
# ---------------------------------------------------------------------------
AROMATIC_COMPOUND_NAMES: Sequence[str] = (
    # 4-F-Phe
    "4-fluorophenylalanine", "para-fluorophenylalanine",
    "p-fluorophenylalanine", "4-FPhe",
    # 4-Cl-Phe
    "4-chlorophenylalanine", "para-chlorophenylalanine",
    "p-chlorophenylalanine",
    # 3-Cl-Phe
    "3-chlorophenylalanine", "meta-chlorophenylalanine",
    # 2-F-Phe
    "2-fluorophenylalanine", "ortho-fluorophenylalanine",
    # 3-Cl,4-F-Phe
    "3-chloro-4-fluorophenylalanine",
    # 4-CN-Phe
    "4-cyanophenylalanine", "para-cyanophenylalanine",
    # 4-Me-Phe
    "4-methylphenylalanine", "para-methylphenylalanine",
    # 3-Me-Phe
    "3-methylphenylalanine",
    # 4-CF3-Phe
    "4-trifluoromethylphenylalanine",
    "para-trifluoromethylphenylalanine",
    # Bip — biphenylalanine
    "biphenylalanine", "4-biphenylalanine",
    "4-phenylphenylalanine", "Bip",
    "3-(4-biphenylyl)alanine",
    # 2-Nal — 2-naphthylalanine
    "2-naphthylalanine", "L-2-naphthylalanine",
    "beta-naphthylalanine", "3-(2-naphthyl)alanine",
    "3-(2-naphthyl)-L-alanine",
    # 1-Nal — 1-naphthylalanine
    "1-naphthylalanine", "L-1-naphthylalanine",
    "alpha-naphthylalanine", "3-(1-naphthyl)alanine",
    # 2-Pal / 3-Pal / 4-Pal — pyridylalanine
    "2-pyridylalanine", "3-pyridylalanine", "4-pyridylalanine",
    "pyridylalanine",
    "3-(2-pyridyl)alanine", "3-(3-pyridyl)alanine", "3-(4-pyridyl)alanine",
    # 4-OMe-Phe — 4-methoxyphenylalanine
    "4-methoxyphenylalanine", "O-methyltyrosine",
    "para-methoxyphenylalanine",
    # 4-NH2-Phe — 4-aminophenylalanine
    "4-aminophenylalanine", "para-aminophenylalanine",
    # 3-F-Tyr — 3-fluorotyrosine
    "3-fluorotyrosine", "3-fluoro-L-tyrosine",
    # Bta — benzothienylalanine
    "benzothienylalanine", "beta-benzothienylalanine",
    "3-(2-benzothienyl)alanine", "3-(3-benzothienyl)alanine",
    # 3-Thi — thienylalanine
    "thienylalanine", "3-thienylalanine",
    "3-(thiophen-3-yl)alanine", "beta-thienylalanine",
    # h-Phe — homophenylalanine
    "homophenylalanine", "L-homophenylalanine",
    "2-amino-4-phenylbutanoic acid",
    # Tic — 1,2,3,4-tetrahydroisoquinoline-3-carboxylic acid
    "1,2,3,4-tetrahydroisoquinoline-3-carboxylic acid",
    "Tic amino acid", "L-Tic",
    "tetrahydroisoquinolinecarboxylic acid",
    # Igl — indanylglycine
    "indanylglycine",
    # 3-Me-His — 3-methylhistidine
    "3-methylhistidine", "tau-methylhistidine",
    # 4-CONH2-Phg
    "4-carbamoylphenylglycine",
    # Phe[4-(2-aminoethoxy)]
    "4-(2-aminoethoxy)phenylalanine",
    # 33DPA — 3,3-diphenylalanine
    "3,3-diphenylalanine", "alpha,alpha-diphenylalanine",
)

# ---------------------------------------------------------------------------
# Pro / 脯氨酸类（环状亚氨基酸）
# ---------------------------------------------------------------------------
PRO_COMPOUND_NAMES: Sequence[str] = (
    # trans-HyP / Hyp — 4-hydroxyproline
    "4-hydroxyproline", "trans-4-hydroxyproline",
    "L-4-hydroxyproline", "(2S,4R)-4-hydroxyproline",
    "trans-L-hydroxyproline", "hydroxyproline",
    # Hyp(3-OH) — 3-hydroxyproline
    "3-hydroxyproline", "cis-3-hydroxyproline",
    # Pro(4-F2) — 4,4-difluoroproline
    "4,4-difluoroproline", "difluoroproline",
    # Hyp(Et), Hyp(Bzl) — O-substituted hydroxyproline
    "O-ethylhydroxyproline", "O-benzylhydroxyproline",
    # Pro(4-keto) — 4-oxoproline / 4-ketoproline
    "4-oxoproline", "4-ketoproline",
    # Pro(4-Ph) — 4-phenylproline
    "4-phenylproline", "trans-4-phenyl-L-proline",
    # Pro(4R-Et), Pro(4R-nPr) — 4-substituted prolines
    "4-ethylproline", "(4R)-4-ethyl-L-proline",
    "4-propylproline", "(4R)-4-propyl-L-proline",
    # hPro — homoproline
    "homoproline", "L-homoproline", "pipecolic acid",
    "pipecolinic acid", "L-pipecolic acid",
    "2-piperidinecarboxylic acid",
    # h-Pic — homopipecolic acid
    "homopipecolic acid", "nipecotic acid homologue",
    # Aze — azetidine-2-carboxylic acid
    "azetidine-2-carboxylic acid", "L-azetidine-2-carboxylic acid",
    "azetidine-2-carboxylic acid", "Aze",
    # Pip/Pic — pipecolic acid (already covered in hPro above)
    "pipecolic acid", "pipecolinic acid",
    # Pic(4-Oxo) — 4-oxopipecolic acid
    "4-oxopipecolic acid", "4-ketopipecolic acid",
    # Thz / ThioPro — thiaproline
    "thiazolidine-4-carboxylic acid", "thiaproline", "thioproline",
    "L-thiazolidine-4-carboxylic acid",
    # Aze(3-Me), Aze(3-Me2) — methylazetidine
    "3-methylazetidine-2-carboxylic acid",
    "3,3-dimethylazetidine-2-carboxylic acid",
    # Tic — tetrahydroisoquinoline (also aromatic)
    "1,2,3,4-tetrahydroisoquinoline-3-carboxylic acid",
    # Idc — indoline-2-carboxylic acid
    "indoline-2-carboxylic acid", "L-indoline-2-carboxylic acid",
    # Oic — octahydroindole-2-carboxylic acid
    "octahydroindole-2-carboxylic acid",
    "(3aS,7aS)-octahydroindole-2-carboxylic acid",
    "Oic", "perhydroindole-2-carboxylic acid",
    # Mor — morpholine-3-carboxylic acid
    "morpholine-3-carboxylic acid", "morpholino amino acid",
)

# ---------------------------------------------------------------------------
# Trp / 吲哚类衍生
# ---------------------------------------------------------------------------
TRP_COMPOUND_NAMES: Sequence[str] = (
    # 5-F-Trp
    "5-fluorotryptophan", "5-fluoro-L-tryptophan",
    # 7-F-Trp
    "7-fluorotryptophan", "7-fluoro-L-tryptophan",
    # Me-Trp, 7-Me-Trp, 2-Me-Trp
    "5-methyltryptophan", "7-methyltryptophan", "2-methyltryptophan",
    "methyltryptophan", "N-methyltryptophan",
    # beta-Me-Trp
    "beta-methyltryptophan", "alpha-methyltryptophan",
    # 5-OH-Trp
    "5-hydroxytryptophan", "5-HTP", "oxitriptan",
    "5-hydroxy-L-tryptophan",
    # AzaTrp
    "azatryptophan", "aza-tryptophan",
    # Dht — dehydrotryptophan
    "dehydrotryptophan", "2,3-dehydrotryptophan",
    "dihydrotryptophan",
    # Bta — benzothienylalanine (shared with aromatic)
    "benzothienylalanine",
)

# ---------------------------------------------------------------------------
# Ser / Thr / 含羟基衍生
# ---------------------------------------------------------------------------
SER_COMPOUND_NAMES: Sequence[str] = (
    # Ser(Me) — O-methylserine
    "O-methylserine", "O-methyl-L-serine",
    # Ser(Ac) — O-acetylserine
    "O-acetylserine", "O-acetyl-L-serine",
    # Ser(Bzl) — O-benzylserine
    "O-benzylserine", "O-benzyl-L-serine",
    # Ser(iPen) — O-isopentylserine
    "O-isopentylserine",
    # Ser(nPr) — O-propylserine
    "O-propylserine", "O-n-propylserine",
    # Ser(cPr) — O-cyclopropylserine
    "O-cyclopropylserine",
    # Ser(cBu) — O-cyclobutylserine
    "O-cyclobutylserine",
    # Thr(Me) — O-methylthreonine
    "O-methylthreonine",
    # hSer(Ph-4-Cl)
    "O-(4-chlorobenzyl)homoserine",
    # Ala(B(OH)2) — borono-alanine
    "borono-alanine", "boronoalanine",
    "2-amino-3-boronopropanoic acid",
    "L-boronoalanine",
    # Val(3-OH) — beta-hydroxyvaline
    "3-hydroxyvaline", "beta-hydroxyvaline",
    "3-hydroxy-L-valine",
    # Ser(TFE) — trifluoroethylserine
    "O-(2,2,2-trifluoroethyl)serine", "trifluoroethylserine",
)

# ---------------------------------------------------------------------------
# N-取代甘氨酸（肽拟肽骨架）
# ---------------------------------------------------------------------------
GLY_PEPTOID_NAMES: Sequence[str] = (
    "N-cyclopropylglycine",
    "N-isopropylglycine",
    "N-propylglycine",
    "N-(4-chlorobenzyl)glycine",
    "N-thienylglycine", "N-(thiophen-3-yl)glycine",
    "N-(2-pyridinylmethyl)glycine",
    "N-(4-pyridinylmethyl)glycine",
    "N-benzylglycine", "sarcosine",
    "N-(2-phenylethyl)glycine",
    "N-cyclopentylglycine",
)

# ---------------------------------------------------------------------------
# Asp/Asn/Gln 酰胺衍生类
# ---------------------------------------------------------------------------
ASP_ASN_COMPOUND_NAMES: Sequence[str] = (
    "aspartyl piperidine amide",
    "aspartyl morpholine amide",
    "asparagine piperidine amide",
    "asparagine morpholine amide",
    "glutamine piperidine amide",
    "glutamine morpholine amide",
    "asparagine azetidine amide",
    "glutamine azetidine amide",
    "difluoropiperidyl asparagine",
    "difluoropiperidyl glutamine",
    "aspartamide amino acid",
)

# ---------------------------------------------------------------------------
# β-氨基酸类
# ---------------------------------------------------------------------------
BETA_AA_NAMES: Sequence[str] = (
    # β-Ala — beta-alanine（本体，非蛋白生成氨基酸）
    "beta-alanine", "3-aminopropionic acid", "beta-aminopropionic acid",
    # 3-Abu — beta-aminobutyric acid
    "3-aminobutyric acid", "beta-aminobutyric acid",
    "3-aminobutanoic acid",
    # 3-(CF3)-bAla — trifluoromethyl beta-alanine
    "3-(trifluoromethyl)-beta-alanine",
    "3-trifluoromethyl-beta-alanine",
    # 2-ACHxC — trans-2-aminocyclohexane carboxylic acid
    "2-aminocyclohexanecarboxylic acid",
    "trans-2-aminocyclohexane carboxylic acid",
    "trans-ACHC",
    # 2-ACPnC — 2-aminocyclopentane carboxylic acid
    "2-aminocyclopentanecarboxylic acid",
    "trans-2-aminocyclopentane carboxylic acid",
    # R-AMPA — alpha-methylphenylalanine
    "alpha-methylphenylalanine", "R-alpha-methylphenylalanine",
    "2-amino-2-methyl-3-phenylpropionic acid",
    # h-Hph — homohomophenylalanine
    "homohomophenylalanine", "2-amino-5-phenylpentanoic acid",
    # general
    "beta-amino acid synthesis",
)

# ---------------------------------------------------------------------------
# GCE 专属化合物（来自 iNCLusive 高频 ncAA，遗传密码扩展中最常见的化合物）
# 这些化合物的合成/应用/代谢通路均在 GCE 轨道中检索
# ---------------------------------------------------------------------------
GCE_COMPOUND_NAMES: Sequence[str] = (
    # ── 苯丙氨酸衍生物（Phe-based, click chemistry / photo-crosslinking）──
    # pAzF / AzF — p-azidophenylalanine（最常见 GCE ncAA）
    "p-azidophenylalanine", "para-azidophenylalanine",
    "4-azidophenylalanine", "4-azido-L-phenylalanine",
    "p-azido-L-phenylalanine", "pAzF", "AzF",
    # pBpa / Bpa — p-benzoylphenylalanine（光交联探针）
    "p-benzoylphenylalanine", "para-benzoylphenylalanine",
    "4-benzoylphenylalanine", "4-benzoyl-L-phenylalanine",
    "p-benzoyl-L-phenylalanine", "pBpa", "Bpa",
    # pAcF / pAcPhe / AcF — p-acetylphenylalanine（生物正交 oxime/hydrazone）
    "p-acetylphenylalanine", "para-acetylphenylalanine",
    "4-acetylphenylalanine", "4-acetyl-L-phenylalanine",
    "p-acetyl-L-phenylalanine", "pAcF", "pAcPhe", "AcF",
    # pNO2F / pNBF — p-nitrophenylalanine
    "p-nitrophenylalanine", "4-nitrophenylalanine",
    "para-nitrophenylalanine", "4-nitro-L-phenylalanine",
    # pIF — p-iodophenylalanine
    "p-iodophenylalanine", "4-iodophenylalanine",
    "para-iodophenylalanine",
    # pBrF — p-bromophenylalanine
    "p-bromophenylalanine", "4-bromophenylalanine",
    # pCNF — p-cyanophenylalanine (same as 4-CN-Phe above, add abbreviation)
    "pCNF",
    # MeaF — meta-aminophenylalanine
    "meta-aminophenylalanine", "3-aminophenylalanine",
    # AzAla — beta-azidoalanine
    "azidoalanine", "beta-azidoalanine", "2-amino-3-azidopropanoic acid",
    # ── 酪氨酸衍生物（Tyr-based）──
    # nitroTyr — 3-nitrotyrosine（广泛用于 GCE 和蛋白质修饰研究）
    "3-nitrotyrosine", "nitrotyrosine", "3-nitro-L-tyrosine",
    "nitroTyr", "L-3-nitrotyrosine",
    # DOPA — L-DOPA / 3,4-dihydroxyphenylalanine
    "3,4-dihydroxyphenylalanine", "L-DOPA", "levodopa",
    "3,4-dihydroxy-L-phenylalanine", "DOPA",
    # ONBY — o-nitrobenzyl tyrosine（光脱保护 ncAA）
    "o-nitrobenzyltyrosine", "ortho-nitrobenzyltyrosine",
    "2-nitrobenzyltyrosine", "ONBY",
    # BetY / sTyr — sulfotyrosine
    "sulfotyrosine", "O-sulfo-L-tyrosine", "sTyr",
    # OpgY — O-phosphoglycol tyrosine
    "phosphoglycol tyrosine",
    # ── 赖氨酸衍生物（Lys-based, GCE epigenetics / bioorthogonal）──
    # BocK — Nε-Boc-L-lysine
    "N-epsilon-Boc-lysine", "Nepsilon-tert-butoxycarbonyllysine",
    "Boc-lysine", "BocK", "Nε-Boc-lysine",
    # AllocK / AlocK — Nε-allyloxycarbonyl-L-lysine
    "N6-allyloxycarbonyllysine", "N-epsilon-allyloxycarbonyllysine",
    "N6-[(allyloxy)carbonyl]-L-lysine", "AllocK", "AlocK",
    # PrK / ProK — propargyllysine / Nε-propargyloxycarbonyl-lysine
    "propargyllysine", "N-epsilon-propargyllysine",
    "Nε-propargyloxycarbonyl-lysine", "PrK", "ProK",
    # AcK — Nε-acetyllysine（also in LYS list, add abbreviation）
    "AcK",
    # BCNK — bicyclononyne lysine
    "bicyclononyne lysine", "Nε-bicyclononyneoxycarbonyl-lysine",
    "BCNK",
    # TCO*K — trans-cyclooctene lysine
    "trans-cyclooctene lysine", "TCO*K", "TCO-lysine",
    # AzK — azide lysine
    "azide lysine", "Nε-azidolysine", "AzK",
    # CrK — crotonyllysine
    "crotonyllysine", "Nε-crotonyl-lysine", "CrK",
    # ThrK — threonyl lysine
    "threonyl lysine", "ThrK",
    # VtK — vinyl thioether lysine
    "VtK",
    # NAEK — Nε-azidoethyloxycarbonyllysine
    "NAEK",
    # AbK — AbK
    "AbK",
    # DiZPK
    "diazidopropionyl lysine", "DiZPK",
    # TMSK
    "TMSK",
    # ── 丝氨酸/苏氨酸衍生物（phosphorylation signals）──
    # Sep — phosphoserine
    "phosphoserine", "O-phosphoserine", "L-phosphoserine",
    "O-phospho-L-serine", "Sep",
    # pSer — phosphoserine abbreviation
    "pSer",
    # ── 色氨酸衍生物（Trp-based）──
    # Acd — acridonylalanine
    "acridonylalanine", "Acd",
    # ── 荧光/功能性 ncAA ──
    # Cou — coumaryl amino acid
    "coumaryl amino acid", "coumarinaianine", "Cou",
    # ANL — azidonorleucine（甲硫氨酸类似物，广泛用于 BONCAT）
    "azidonorleucine", "L-azidonorleucine",
    "2-amino-6-azidohexanoic acid", "ANL",
    # Azi — azidophenylalanine variants
    "Azi",
    # ── 更多 Phe 衍生物（bioorthogonal，来自 iNCLusive 高频数据）──
    # pCNF — p-cyanophenylalanine
    "p-cyanophenylalanine", "4-cyanophenylalanine", "para-cyanophenylalanine",
    "4-cyano-L-phenylalanine", "pCNF",
    # pBoF — p-boronophenylalanine
    "p-boronophenylalanine", "4-boronophenylalanine", "para-boronophenylalanine",
    "4-borono-L-phenylalanine", "pBoF",
    # pPrOF — p-propargyloxyphenylalanine
    "p-propargyloxyphenylalanine", "4-propargyloxyphenylalanine",
    "O-propargyl-tyrosine", "O-propargyl-L-tyrosine", "pPrOF", "OPgY",
    # ClAcF — p-chloroacetamidophenylalanine
    "p-chloroacetamidophenylalanine", "4-chloroacetamidophenylalanine",
    # 3-Amino-Tyr / 3-AminoTyr
    "3-aminotyrosine", "3-amino-L-tyrosine",
    # o-NB-Glu — o-nitrobenzyl-glutamic acid (photocaged)
    "o-nitrobenzyl-L-glutamic acid", "O-nitrobenzyl glutamic acid",
    # ── 更多 Lys 衍生物（iNCLusive 高频）──
    # PCK — pyrroline-carboxy-lysine (pyrrolysyl analog)
    "pyrroline-carboxy-lysine", "PCK", "pyrrolysine analog",
    # coumarin-Lys / aminocoumarin lysine
    "aminocoumarin lysine", "coumarin lysine",
    "L-(7-hydroxycoumarin-4-yl)ethylglycine", "7-hydroxycoumarin alanine",
    # NAEK variants / azidoethyloxycarbonyllysine
    "Nε-[(2-azideoethyloxy)carbonyl]-L-lysine", "azidoethyloxycarbonyllysine",
    # DiZ-diazirine-Lys
    "Nε-[(2-(3-methyl-3H-diazirin-3-yl)ethoxy)carbonyl]-L-lysine",
    "diazirine lysine", "photo-crosslinking lysine",
    # Axial TCO-Lys
    "axial trans-cyclooct-2-ene-L-lysine", "axial TCO-lysine",
    # VtK — vinyl thioether lysine
    "N6-((2-(vinylthio)ethoxy)carbonyl)-L-lysine",
    # alkynyl-Lys
    "N-propargyl-L-lysine", "propargyl-L-lysine",
    "alkynyllysine", "alkynyl lysine",
    # ── 更多 Tyr 衍生物（iNCLusive 高频）──
    # ortho-NB-Tyr
    "ortho-nitrobenzyl tyrosine", "O-nitrobenzyl-L-tyrosine",
    # o-bromo-Tyr
    "O-(2-bromoethyl)-L-tyrosine",
    # ── Cys 衍生物 ──
    # S-allyl-Cys
    "S-allyl-L-cysteine", "S-allylcysteine",
    # S-propargyl-Cys
    "S-propargyl-L-cysteine",
    # ── 甲硫氨酸类似物 ──
    "azidohomoalanine", "L-azidohomoalanine", "AHA",
    "homopropargylglycine", "L-homopropargylglycine", "HPG",
    # ── 通用 GCE 相关化合物类别 ──
    "pyrrolysine", "Pyl", "L-pyrrolysine",
    "(2R)-2-amino-3-(((2S,3R)-3-methyl-2-[(1-oxopyrrolidin-2-ylidene)amino]butanoyl)oxy)propanoic acid",
    "selenocysteine", "selenocysteine amino acid",
)

# ===========================================================================
# 第二部分：合并全部化合物名 → 用于 CrossRef
# ===========================================================================

ALL_PDF_COMPOUND_NAMES: Sequence[str] = (
    *LYS_COMPOUND_NAMES,
    *ARG_COMPOUND_NAMES,
    *ALIPHATIC_COMPOUND_NAMES,
    *CYCLIC_COMPOUND_NAMES,
    *AROMATIC_COMPOUND_NAMES,
    *PRO_COMPOUND_NAMES,
    *TRP_COMPOUND_NAMES,
    *SER_COMPOUND_NAMES,
    *GLY_PEPTOID_NAMES,
    *ASP_ASN_COMPOUND_NAMES,
    *BETA_AA_NAMES,
    *GCE_COMPOUND_NAMES,
)

# ===========================================================================
# 第三部分：PubMed 检索块（按化合物组）
# ===========================================================================

def _pm_block(names: Sequence[str]) -> str:
    """将化合物名列表拼接为 PubMed Title/Abstract OR 块。"""
    return " OR ".join(f"{n}[Title/Abstract]" for n in names)


_PM_LYS      = _pm_block(LYS_COMPOUND_NAMES)
_PM_ARG      = _pm_block(ARG_COMPOUND_NAMES)
_PM_ALIPH    = _pm_block(ALIPHATIC_COMPOUND_NAMES)
_PM_CYCLIC   = _pm_block(CYCLIC_COMPOUND_NAMES)
_PM_AROM     = _pm_block(AROMATIC_COMPOUND_NAMES)
_PM_PRO      = _pm_block(PRO_COMPOUND_NAMES)
_PM_TRP      = _pm_block(TRP_COMPOUND_NAMES)
_PM_SER      = _pm_block(SER_COMPOUND_NAMES)
_PM_BETA     = _pm_block(BETA_AA_NAMES)

# 合并高频化合物块（用于 AI 筛选 prompt / 宽泛检索 fallback）
_PM_ALL_COMPOUNDS = _pm_block((
    "norvaline", "norleucine", "tert-leucine", "Aib",
    "alpha-aminoisobutyric acid", "cyclohexylalanine",
    "4-fluorophenylalanine", "2-naphthylalanine", "1-naphthylalanine",
    "biphenylalanine", "4-hydroxyproline", "hydroxyproline",
    "pipecolic acid", "citrulline", "homoarginine", "ornithine",
    "5-hydroxytryptophan", "5-fluorotryptophan",
    "azetidine-2-carboxylic acid", "4,4-difluoroproline",
    "diaminopropionic acid", "diaminobutyric acid",
    "homophenylalanine", "propargylglycine", "allylglycine",
    "phenylglycine", "cyclohexylglycine", "homoserine",
))

# ===========================================================================
# 第四部分：各轨道 PubMed 检索式
# ===========================================================================

# --- 工艺信号词 ---
_CHEM_SIG   = ("chemical synthesis[Title/Abstract] OR asymmetric synthesis[Title/Abstract] "
               "OR enantioselective synthesis[Title/Abstract] OR stereoselective synthesis[Title/Abstract] "
               "OR total synthesis[Title/Abstract] OR organocatalysis[Title/Abstract]")
_ENZ_SIG    = ("enzymatic synthesis[Title/Abstract] OR biocatalysis[Title/Abstract] "
               "OR transaminase[Title/Abstract] OR aminotransferase[Title/Abstract] "
               "OR reductive amination[Title/Abstract] OR enzyme cascade[Title/Abstract]")
_FERM_SIG   = ("fermentation[Title/Abstract] OR microbial production[Title/Abstract] "
               "OR metabolic engineering[Title/Abstract] OR cell factory[Title/Abstract] "
               "OR whole-cell biocatalysis[Title/Abstract]")
_PATH_SIG   = ("biosynthetic pathway[Title/Abstract] OR metabolic pathway[Title/Abstract] "
               "OR biosynthesis[Title/Abstract] OR de novo biosynthesis[Title/Abstract] "
               "OR pathway engineering[Title/Abstract]")
_HYB_SIG    = ("chemoenzymatic[Title/Abstract] OR chemoenzymatic synthesis[Title/Abstract] "
               "OR sequential biocatalysis[Title/Abstract] OR one-pot biocatalysis[Title/Abstract] "
               "OR enzyme-chemical cascade[Title/Abstract] "
               "OR chemoenzymatic route[Title/Abstract] OR chemo-enzymatic[Title/Abstract]")

# 广义 NNAA 术语——用于捕捉不点名具体化合物的方法论/综合性文献
# 注意：这些检索式召回范围更宽，依赖 AI 筛选把关
_NNAA_GENERIC = (
    "non-natural amino acid[Title/Abstract] OR unnatural amino acid[Title/Abstract] "
    "OR noncanonical amino acid[Title/Abstract] OR non-canonical amino acid[Title/Abstract] "
    "OR ncAA[Title/Abstract]"
)

# --- GCE 遗传密码扩展专属信号词 ---
# aaRS 系统名（最高特异性，几乎仅出现在 GCE 论文中）
_GCE_AARS = (
    "PylRS[Title/Abstract] OR MjTyrRS[Title/Abstract] OR EcTyrRS[Title/Abstract] "
    "OR MmPylRS[Title/Abstract] OR MbPylRS[Title/Abstract] "
    "OR pyrrolysyl-tRNA synthetase[Title/Abstract] "
    "OR orthogonal aminoacyl-tRNA synthetase[Title/Abstract]"
)
# GCE 方法学词（通用性稍宽，但与 ncAA 组合精确）
_GCE_METHOD = (
    "genetic code expansion[Title/Abstract] OR expanded genetic code[Title/Abstract] "
    "OR amber suppression[Title/Abstract] OR stop codon suppression[Title/Abstract] "
    "OR orthogonal tRNA[Title/Abstract] OR unnatural amino acid incorporation[Title/Abstract] "
    "OR noncanonical amino acid incorporation[Title/Abstract] "
    "OR non-canonical amino acid incorporation[Title/Abstract]"
)
# GCE 高频化合物块
_PM_GCE_COMPOUNDS = _pm_block((
    "p-azidophenylalanine", "p-benzoylphenylalanine", "p-acetylphenylalanine",
    "4-azidophenylalanine", "4-benzoylphenylalanine", "4-acetylphenylalanine",
    "phosphoserine", "3-nitrotyrosine", "pyrrolysine",
    "azidonorleucine", "propargyllysine", "trans-cyclooctene lysine",
    "DOPA", "3,4-dihydroxyphenylalanine",
))

# --- 代谢通路（pathway） ---
PUBMED_PATHWAY_QUERIES: Sequence[str] = (
    f"({_PM_ALIPH}) AND ({_PATH_SIG})",
    f"({_PM_CYCLIC}) AND ({_PATH_SIG})",
    f"({_PM_ARG}) AND ({_PATH_SIG})",
    f"({_PM_LYS}) AND ({_PATH_SIG})",
    f"({_PM_AROM}) AND ({_PATH_SIG})",
    f"({_PM_PRO}) AND ({_PATH_SIG})",
    f"({_PM_TRP}) AND ({_PATH_SIG})",
    f"({_PM_SER}) AND ({_PATH_SIG})",
    # 高频化合物宽泛通路检索
    f"({_PM_ALL_COMPOUNDS}) AND ({_PATH_SIG})",
    # 广义 NNAA + 通路（捕捉不列具体化合物名的代谢通路文献）+ 限定英文
    f"({_NNAA_GENERIC}) AND ({_PATH_SIG}) AND English[Language]",
)

# --- 酶法合成（enzymatic） ---
PUBMED_ENZYMATIC_QUERIES: Sequence[str] = (
    f"({_PM_ALIPH}) AND ({_ENZ_SIG})",
    f"({_PM_CYCLIC}) AND ({_ENZ_SIG})",
    f"({_PM_AROM}) AND ({_ENZ_SIG})",
    f"({_PM_PRO}) AND ({_ENZ_SIG})",
    f"({_PM_ARG}) AND ({_ENZ_SIG})",
    f"({_PM_LYS}) AND ({_ENZ_SIG})",
    f"({_PM_TRP}) AND ({_ENZ_SIG})",
    f"({_PM_BETA}) AND ({_ENZ_SIG})",
    # 高频化合物 + 酶法
    f"({_PM_ALL_COMPOUNDS}) AND ({_ENZ_SIG})",
    # 广义 NNAA + 酶法（捕捉转氨酶、酶级联等方法论文章）
    f"({_NNAA_GENERIC}) AND ({_ENZ_SIG})",
)

# --- 生物发酵（fermentation） ---
PUBMED_FERMENTATION_QUERIES: Sequence[str] = (
    f"({_PM_ALIPH}) AND ({_FERM_SIG})",
    f"({_PM_ARG}) AND ({_FERM_SIG})",
    f"({_PM_LYS}) AND ({_FERM_SIG})",
    f"({_PM_AROM}) AND ({_FERM_SIG})",
    f"({_PM_PRO}) AND ({_FERM_SIG})",
    f"({_PM_TRP}) AND ({_FERM_SIG})",
    f"({_PM_SER}) AND ({_FERM_SIG})",
    # 高频化合物 + 发酵
    f"({_PM_ALL_COMPOUNDS}) AND ({_FERM_SIG})",
    # 广义 NNAA + 发酵（捕捉细胞工厂/代谢工程生产 NNAA 的综合性文献）
    f"({_NNAA_GENERIC}) AND ({_FERM_SIG})",
)

# --- 化学合成（chemical） ---
PUBMED_CHEMICAL_QUERIES: Sequence[str] = (
    # 各结构类 + 化学合成
    f"({_PM_ALIPH}) AND ({_CHEM_SIG})",
    f"({_PM_CYCLIC}) AND ({_CHEM_SIG})",
    f"({_PM_AROM}) AND ({_CHEM_SIG})",
    f"({_PM_PRO}) AND ({_CHEM_SIG})",
    f"({_PM_ARG}) AND ({_CHEM_SIG})",
    f"({_PM_LYS}) AND ({_CHEM_SIG})",
    f"({_PM_TRP}) AND ({_CHEM_SIG})",
    f"({_PM_SER}) AND ({_CHEM_SIG})",
    f"({_PM_BETA}) AND ({_CHEM_SIG})",
    # 高频化合物 + 化学
    f"({_PM_ALL_COMPOUNDS}) AND ({_CHEM_SIG})",
    # 广义 NNAA + 化学合成
    f"({_NNAA_GENERIC}) AND ({_CHEM_SIG})",
    # 芳香类特定组合（文献量大，单独列出）
    "4-fluorophenylalanine[Title/Abstract] AND (synthesis[Title/Abstract] OR enantioselective[Title/Abstract])",
    "2-naphthylalanine[Title/Abstract] AND (synthesis[Title/Abstract] OR preparation[Title/Abstract])",
    "biphenylalanine[Title/Abstract] AND (synthesis[Title/Abstract])",
    "cyclohexylalanine[Title/Abstract] AND (synthesis[Title/Abstract] OR asymmetric[Title/Abstract])",
    "tert-leucine[Title/Abstract] AND (synthesis[Title/Abstract] OR resolution[Title/Abstract])",
    "norvaline[Title/Abstract] AND (synthesis[Title/Abstract] OR preparation[Title/Abstract])",
    "4-hydroxyproline[Title/Abstract] AND (synthesis[Title/Abstract])",
    "pipecolic acid[Title/Abstract] AND (synthesis[Title/Abstract])",
)

# --- 酶-化学联用（hybrid） ---
PUBMED_HYBRID_QUERIES: Sequence[str] = (
    f"({_PM_ALIPH}) AND ({_HYB_SIG})",
    f"({_PM_AROM}) AND ({_HYB_SIG})",
    f"({_PM_CYCLIC}) AND ({_HYB_SIG})",
    f"({_PM_PRO}) AND ({_HYB_SIG})",
    f"({_PM_ARG}) AND ({_HYB_SIG})",
    f"({_PM_ALL_COMPOUNDS}) AND ({_HYB_SIG})",
    # 广义 NNAA + 联用合成
    f"({_NNAA_GENERIC}) AND ({_HYB_SIG})",
)

# --- GCE 遗传密码扩展（gce）---
# 覆盖 aaRS/tRNA 系统介导的 ncAA 蛋白质整合、应用及相关合成研究
PUBMED_GCE_QUERIES: Sequence[str] = (
    # 方法学词（最高精度）
    f"({_GCE_METHOD})",
    # aaRS 系统名（最高特异性）
    f"({_GCE_AARS})",
    # GCE 高频化合物 + 方法学词组合
    f"({_PM_GCE_COMPOUNDS}) AND ({_GCE_METHOD})",
    # GCE 化合物 + aaRS/tRNA 信号
    f"({_PM_GCE_COMPOUNDS}) AND ({_GCE_AARS})",
    # 广义 NNAA + GCE 方法
    f"({_NNAA_GENERIC}) AND ({_GCE_METHOD})",
    # 具体化合物独立检索（高频化合物，补充未用 GCE 术语但明确是 GCE 工作的论文）
    "p-azidophenylalanine[Title/Abstract]",
    "p-benzoylphenylalanine[Title/Abstract]",
    "p-acetylphenylalanine[Title/Abstract]",
    "pyrrolysine[Title/Abstract]",
    "phosphoserine[Title/Abstract] AND (protein incorporation[Title/Abstract] OR aaRS[Title/Abstract] OR amber[Title/Abstract])",
    "3-nitrotyrosine[Title/Abstract] AND (genetic code[Title/Abstract] OR protein[Title/Abstract])",
    "azidonorleucine[Title/Abstract]",
)

# ===========================================================================
# 第五部分：CrossRef 信号词（按轨道）
# ===========================================================================

# 各轨道工艺信号词（CrossRef 短语检索）
_CR_CHEM_METHODS: Sequence[str] = (
    "asymmetric synthesis",
    "enantioselective synthesis",
    "stereoselective synthesis",
    "total synthesis",
    "chemical synthesis",
    "organocatalysis",
    "phase-transfer catalysis",
    "photoredox catalysis",
    "palladium-catalyzed",
    "chiral auxiliary",
    "Strecker synthesis",
    "dynamic kinetic resolution",
    "kinetic resolution",
    "C-H amination",
)

_CR_ENZ_METHODS: Sequence[str] = (
    "enzymatic synthesis",
    "biocatalytic synthesis",
    "biocatalysis",
    "transaminase",
    "aminotransferase",
    "reductive amination",
    "enzyme cascade",
    "multi-enzyme cascade",
    "cell-free protein synthesis",
    "aminoacyl-tRNA synthetase",
    "hydantoinase",
    "carbamoylase",
    "ammonia lyase",
)

_CR_FERM_METHODS: Sequence[str] = (
    "fermentative production",
    "microbial production",
    "microbial synthesis",
    "metabolic engineering amino acid",
    "cell factory amino acid",
    "whole-cell biocatalysis",
    "recombinant strain amino acid",
    "fed-batch fermentation amino acid",
)

_CR_PATH_METHODS: Sequence[str] = (
    "biosynthetic pathway",
    "amino acid biosynthesis",
    "de novo biosynthesis",
    "metabolic pathway engineering",
    "pathway engineering",
    "precursor supply",
    "flux redistribution",
    "transporter engineering",
)

_CR_HYB_METHODS: Sequence[str] = (
    "chemoenzymatic synthesis",
    "chemoenzymatic preparation",
    "sequential chemoenzymatic",
    "one-pot biocatalysis",
    "enzyme chemical cascade",
    "chemoenzymatic route",
)

# GCE 专属 CrossRef 检索词（精度极高）
_CR_GCE_METHODS: Sequence[str] = (
    "genetic code expansion",
    "expanded genetic code",
    "amber suppression",
    "stop codon suppression",
    "orthogonal tRNA synthetase",
    "pyrrolysyl-tRNA synthetase",
    "PylRS",
    "MjTyrRS",
    "unnatural amino acid incorporation",
    "noncanonical amino acid incorporation",
    "site-specific incorporation",
    "orthogonal aaRS",
)

# GCE 高频化合物（CrossRef 短语效果好）
_CR_GCE_COMPOUNDS: Sequence[str] = (
    "p-azidophenylalanine",
    "p-benzoylphenylalanine",
    "p-acetylphenylalanine",
    "4-azidophenylalanine",
    "4-benzoylphenylalanine",
    "phosphoserine",
    "3-nitrotyrosine",
    "pyrrolysine",
    "azidonorleucine",
    "propargyllysine",
    "trans-cyclooctene lysine",
    "3,4-dihydroxyphenylalanine",
)

# 高优先级化合物简名（CrossRef 短语效果好）
_CR_KEY_COMPOUNDS: Sequence[str] = (
    "norvaline",
    "norleucine",
    "tert-leucine",
    "alpha-aminoisobutyric acid",
    "Aib",
    "cyclohexylalanine",
    "4-fluorophenylalanine",
    "fluorophenylalanine",
    "2-naphthylalanine",
    "1-naphthylalanine",
    "naphthylalanine",
    "biphenylalanine",
    "4-hydroxyproline",
    "hydroxyproline",
    "pipecolic acid",
    "citrulline",
    "homoarginine",
    "ornithine",
    "5-hydroxytryptophan",
    "5-fluorotryptophan",
    "fluorotryptophan",
    "azetidine-2-carboxylic acid",
    "difluoroproline",
    "diaminopropionic acid",
    "diaminobutyric acid",
    "homophenylalanine",
    "propargylglycine",
    "allylglycine",
    "phenylglycine",
    "cyclohexylglycine",
    "homoserine",
    "4-pyridylalanine",
    "3-pyridylalanine",
    "2-pyridylalanine",
    "pyridylalanine",
    "4-chlorophenylalanine",
    "chlorophenylalanine",
    "4-cyanophenylalanine",
    "4-methylphenylalanine",
    "4-trifluoromethylphenylalanine",
    "thienylalanine",
    "benzothienylalanine",
    "indanylglycine",
    "cyclopropylalanine",
    "beta-alanine",
    "beta-aminobutyric acid",
    "2-aminocyclohexanecarboxylic acid",
    "homoproline",
    "thiaproline",
    "thiazolidine-4-carboxylic acid",
)

# CrossRef 信号词（per track = 关键化合物 × 工艺词配对）
PATHWAY_SIGNALS_EN: Sequence[str] = (*_CR_KEY_COMPOUNDS, *_CR_PATH_METHODS)
ENZYMATIC_SIGNALS_EN: Sequence[str] = (*_CR_KEY_COMPOUNDS, *_CR_ENZ_METHODS)
FERMENTATION_SIGNALS_EN: Sequence[str] = (*_CR_KEY_COMPOUNDS, *_CR_FERM_METHODS)
CHEMICAL_SIGNALS_EN: Sequence[str] = (*_CR_KEY_COMPOUNDS, *_CR_CHEM_METHODS)
HYBRID_SIGNALS_EN: Sequence[str] = (*_CR_KEY_COMPOUNDS, *_CR_HYB_METHODS)

# ===========================================================================
# 第六部分：汇总索引（供外部模块调用，接口保持兼容）
# ===========================================================================

# compat: NNAA_CORE_EN 保留供 AI 筛选 prompt 引用，限定为 PDF 涉及的化合物类别
NNAA_CORE_EN: Sequence[str] = (
    "non-natural amino acid",
    "unnatural amino acid",
    "noncanonical amino acid",
    "non-proteinogenic amino acid",
    "non-standard amino acid",
)

PUBMED_QUERIES_BY_TRACK: Dict[str, Sequence[str]] = {
    "pathway":      PUBMED_PATHWAY_QUERIES,
    "enzymatic":    PUBMED_ENZYMATIC_QUERIES,
    "fermentation": PUBMED_FERMENTATION_QUERIES,
    "chemical":     PUBMED_CHEMICAL_QUERIES,
    "hybrid":       PUBMED_HYBRID_QUERIES,
    "gce":          PUBMED_GCE_QUERIES,
}

GCE_SIGNALS_EN: Sequence[str] = (*_CR_GCE_COMPOUNDS, *_CR_GCE_METHODS)

CROSSREF_SIGNALS_BY_TRACK: Dict[str, Sequence[str]] = {
    "pathway":      PATHWAY_SIGNALS_EN,
    "enzymatic":    ENZYMATIC_SIGNALS_EN,
    "fermentation": FERMENTATION_SIGNALS_EN,
    "chemical":     CHEMICAL_SIGNALS_EN,
    "hybrid":       HYBRID_SIGNALS_EN,
    "gce":          GCE_SIGNALS_EN,
}


def get_pubmed_queries(track: str) -> List[str]:
    return list(PUBMED_QUERIES_BY_TRACK.get(track, ()))


def _pair_crossref(core_terms: Sequence[str], signals: Sequence[str], max_queries: Optional[int] = None) -> List[str]:
    """化合物名 × 工艺词 笛卡尔积生成 CrossRef 短句。"""
    queries = [f"{core} {signal}".strip() for core, signal in product(core_terms, signals)]
    if max_queries is not None and max_queries > 0:
        return queries[:max_queries]
    return queries


def build_crossref_queries(track: str, max_queries: Optional[int] = None) -> List[str]:
    """
    CrossRef 检索式（title 模式）：
    - 非 GCE 轨道：用 ALL_PDF_COMPOUND_NAMES 全量化合物名直接作为 query.title 检索词，
      同时补充通用 NNAA 术语 + 方法词组合，覆盖化学合成期刊中的方法论文献。
    - GCE 轨道：用 GCE_COMPOUND_NAMES + GCE 方法词。
    所有查询均配合 config 中的 query_param: title 使用，即在标题字段匹配，大幅降低噪声。
    """
    if track == "gce":
        queries: List[str] = list(GCE_COMPOUND_NAMES)
        queries += list(_CR_GCE_METHODS)
        if max_queries and max_queries > 0:
            return queries[:max_queries]
        return queries

    if track not in ("pathway", "enzymatic", "fermentation", "chemical", "hybrid"):
        return []

    # 通用术语 + 各轨道方法词（用于捕捉方法论综述/方法纸）
    method_signals_map = {
        "pathway":      _CR_PATH_METHODS,
        "enzymatic":    _CR_ENZ_METHODS,
        "fermentation": _CR_FERM_METHODS,
        "chemical":     _CR_CHEM_METHODS,
        "hybrid":       _CR_HYB_METHODS,
    }
    method_signals = method_signals_map[track]

    # ① 全量化合物名（548条）直接作为 query.title 词
    queries = list(ALL_PDF_COMPOUND_NAMES)

    # ② NNAA 通用术语 × 方法词（补充方法论论文）
    nnaa_terms = [
        "non-natural amino acid", "unnatural amino acid",
        "noncanonical amino acid", "non-proteinogenic amino acid",
    ]
    for term in nnaa_terms:
        for sig in method_signals:
            queries.append(f"{term} {sig}")

    # ③ hybrid 额外补充联用组合词
    if track == "hybrid":
        queries += [
            "chemoenzymatic synthesis amino acid",
            "one-pot chemoenzymatic amino acid",
            "chemoenzymatic asymmetric amino acid synthesis",
        ]

    if max_queries and max_queries > 0:
        return queries[:max_queries]
    return queries


TRACK_CROSSREF_BUILDERS = {
    track: (lambda t=track: build_crossref_queries(t))
    for track in PUBMED_QUERIES_BY_TRACK
}
