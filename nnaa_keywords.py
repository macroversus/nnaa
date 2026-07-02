"""
非天然氨基酸（NNAA）文献检索关键词库

设计原则：
1. 覆盖范围严格限定于 PDF《Derivative of Natural AA》中列出的全部化合物
2. 每个化合物收录主名 + 所有常见别称/同义词（IUPAC 名、俗名、缩写、立体前缀变体）
3. 按结构母体分组，构建 OR 检索块，再与各轨道工艺词组合
4. 不收录 PDF 未提及的化合物（如 genetic code expansion 等基因编码扩展主题）
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
_HYB_SIG    = ("chemoenzymatic[Title/Abstract] OR one-pot synthesis[Title/Abstract] "
               "OR hybrid synthesis[Title/Abstract] OR sequential biocatalysis[Title/Abstract]")

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
    "metabolic engineering",
    "cell factory",
    "whole-cell biocatalysis",
    "recombinant strain",
    "fed-batch fermentation",
    "Escherichia coli",
    "Corynebacterium glutamicum",
    "Saccharomyces cerevisiae",
)

_CR_PATH_METHODS: Sequence[str] = (
    "biosynthetic pathway",
    "biosynthesis",
    "de novo biosynthesis",
    "metabolic pathway",
    "pathway engineering",
    "precursor supply",
    "flux redistribution",
    "transporter engineering",
)

_CR_HYB_METHODS: Sequence[str] = (
    "chemoenzymatic synthesis",
    "chemoenzymatic",
    "one-pot synthesis",
    "hybrid catalysis",
    "sequential chemoenzymatic",
    "enzyme chemical cascade",
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
}

CROSSREF_SIGNALS_BY_TRACK: Dict[str, Sequence[str]] = {
    "pathway":      PATHWAY_SIGNALS_EN,
    "enzymatic":    ENZYMATIC_SIGNALS_EN,
    "fermentation": FERMENTATION_SIGNALS_EN,
    "chemical":     CHEMICAL_SIGNALS_EN,
    "hybrid":       HYBRID_SIGNALS_EN,
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
    CrossRef 检索式 = 关键化合物名 × 该轨道工艺词
    """
    method_signals_map = {
        "pathway":      _CR_PATH_METHODS,
        "enzymatic":    _CR_ENZ_METHODS,
        "fermentation": _CR_FERM_METHODS,
        "chemical":     _CR_CHEM_METHODS,
        "hybrid":       _CR_HYB_METHODS,
    }
    method_signals = method_signals_map.get(track)
    if not method_signals:
        return []

    # 化合物名 × 工艺词 配对
    queries = _pair_crossref(_CR_KEY_COMPOUNDS, method_signals, max_queries)

    if track == "hybrid":
        extras = [
            "chemoenzymatic synthesis amino acid",
            "one-pot biocatalysis amino acid",
            "sequential enzymatic chemical amino acid",
        ]
        queries.extend(extras)
        if max_queries is not None and max_queries > 0:
            return queries[:max_queries]
    return queries


TRACK_CROSSREF_BUILDERS = {
    track: (lambda t=track: build_crossref_queries(t))
    for track in PUBMED_QUERIES_BY_TRACK
}
