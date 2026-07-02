"""
compound_tagger.py
扫描文献 title / abstract，标注命中的 PDF 化合物名称。

输出字段：
  nnaa_compound   : 逗号分隔的规范化合物名列表，例如 "norvaline, 4-fluorophenylalanine"
                    若未命中任何化合物则为空字符串
"""

from __future__ import annotations

import re
from typing import Dict, List, Sequence, Tuple

# ---------------------------------------------------------------------------
# 别称 → 规范名映射
# 键 = 所有可能在文献中出现的写法（小写）
# 值 = 展示给用户的规范名
# ---------------------------------------------------------------------------
ALIAS_TO_CANONICAL: Dict[str, str] = {
    # --- Lys 衍生 ---
    "d-lysine": "D-lysine",
    "n-epsilon-acetyllysine": "Lys(Ac)",
    "nepsilon-acetyllysine": "Lys(Ac)",
    "n6-acetyllysine": "Lys(Ac)",
    "n,n-dimethyllysine": "Lys(Me2)",
    "n6,n6-dimethyllysine": "Lys(Me2)",
    "dimethyllysine": "Lys(Me2)",
    "trifluoroacetyllysine": "Lys(TFA)",
    "2,3-diaminopropionic acid": "Dap",
    "diaminopropionic acid": "Dap",
    "l-2,3-diaminopropionic acid": "Dap",
    "alpha,beta-diaminopropionic acid": "Dap",
    "2,4-diaminobutyric acid": "Dab",
    "diaminobutyric acid": "Dab",
    "l-2,4-diaminobutyric acid": "Dab",
    "alpha,gamma-diaminobutyric acid": "Dab",
    "ornithine": "Orn",
    "l-ornithine": "Orn",
    "2,5-diaminopentanoic acid": "Orn",

    # --- Arg 衍生 ---
    "homoarginine": "hArg",
    "l-homoarginine": "hArg",
    "n6-amidinolysin": "hArg",
    "n6-(aminoiminomethyl)lysine": "hArg",
    "citrulline": "Cit",
    "l-citrulline": "Cit",
    "2-amino-5-(carbamoylamino)pentanoic acid": "Cit",
    "homocitrulline": "h-Cit",
    "l-homocitrulline": "h-Cit",
    "d-arginine": "D-Arg",
    "norarginine": "norArg",
    "d-norarginine": "D-norArg",
    "d-homoarginine": "D-HomoArg",
    "asymmetric dimethylarginine": "ADMA",
    "adma": "ADMA",
    "ng,ng-dimethylarginine": "ADMA",
    "nomega,nomega-dimethylarginine": "ADMA",
    "monomethylarginine": "MMA",
    "methylarginine": "MMA",
    "ng-monomethylarginine": "MMA",
    "mma": "MMA",
    "ng-methylarginine": "MMA",
    "nitroarginine": "Arg(NO2)",
    "ng-nitro-l-arginine": "Arg(NO2)",
    "symmetric dimethylarginine": "SDMA",
    "sdma": "SDMA",
    "ng,n'g-dimethylarginine": "SDMA",
    "nomega,n'omega-dimethylarginine": "SDMA",
    "cyclopropylarginine": "Arg(c-Pr)",
    "4-guanidinophenylalanine": "Phe(4-guan)",
    "para-guanidinophenylalanine": "Phe(4-guan)",
    "4-guanidinophenylglycine": "Phg(4-guan)",
    "guanidinopropionic acid": "Agp",
    "beta-guanidinopropionic acid": "Agp",
    "guanidinobutyric acid": "Agb",
    "gamma-guanidinobutyric acid": "Agb",

    # --- 脂肪链 ---
    "norvaline": "Nva",
    "l-norvaline": "Nva",
    "2-aminopentanoic acid": "Nva",
    "alpha-aminovaleric acid": "Nva",
    "alpha-aminopentanoic acid": "Nva",
    "norleucine": "Nle",
    "l-norleucine": "Nle",
    "2-aminohexanoic acid": "Nle",
    "alpha-aminocaproic acid": "Nle",
    "tert-leucine": "Tle",
    "l-tert-leucine": "Tle",
    "tert-butylglycine": "Tle",
    "(2s)-2-amino-3,3-dimethylbutanoic acid": "Tle",
    "3,3-dimethyl-l-alanine": "Tle",
    "beta,beta-dimethylnorvaline": "Tle",
    "aib": "Aib",
    "alpha-aminoisobutyric acid": "Aib",
    "2-aminoisobutyric acid": "Aib",
    "alpha-methylalanine": "Aib",
    "alpha-methyl alanine": "Aib",
    "2-amino-2-methylpropanoic acid": "Aib",
    "2-amino-2-methylpropionic acid": "Aib",
    "2-methylalanine": "Aib",
    "2-aminobutyric acid": "Abu",
    "alpha-aminobutyric acid": "Abu",
    "l-2-aminobutyric acid": "Abu",
    "homoleucine": "hLeu",
    "l-homoleucine": "hLeu",
    "neopentylglycine": "NptGly",
    "2-aminoheptanoic acid": "Ahp",
    "alpha-aminoheptanoic acid": "Ahp",
    "2-aminooctanoic acid": "Aoc",
    "methionine sulfone": "Met(O2)",
    "l-methionine sulfone": "Met(O2)",
    "methionine-s,s-dioxide": "Met(O2)",
    "beta-cyanoalanine": "Ala(CN)",
    "3-cyanoalanine": "Ala(CN)",
    "allylglycine": "Algly",
    "l-allylglycine": "Algly",
    "2-amino-4-pentenoic acid": "Algly",
    "propargylglycine": "Pra",
    "l-propargylglycine": "Pra",
    "2-amino-4-pentynoic acid": "Pra",
    "methallylglycine": "MethaGly",
    "3-ethylnorvaline": "Nva(3-Et)",
    "6-hydroxynorleucine": "Nle(6-OH)",
    "3-hydroxyvaline": "Val(3-OH)",
    "beta-hydroxyvaline": "Val(3-OH)",
    "homoserine": "hSer",
    "l-homoserine": "hSer",
    "2-amino-4-hydroxybutyric acid": "hSer",
    "o-methylhomoserine": "hSer(Me)",
    "o-ethylhomoserine": "hSer(Et)",

    # --- 脂环 ---
    "cyclohexylalanine": "Cha",
    "beta-cyclohexylalanine": "Cha",
    "l-cyclohexylalanine": "Cha",
    "3-cyclohexylalanine": "Cha",
    "2-amino-3-cyclohexylpropionic acid": "Cha",
    "homocyclohexylalanine": "h-Cha",
    "4-cyclohexyl-2-aminobutanoic acid": "h-Cha",
    "cyclopropylalanine": "cPrA",
    "beta-cyclopropylalanine": "cPrA",
    "cyclobutylalanine": "Ala(cBu)",
    "cyclopentylalanine": "cPenA",
    "phenylglycine": "Phg",
    "l-phenylglycine": "Phg",
    "d-phenylglycine": "Phg",
    "alpha-phenylglycine": "Phg",
    "aminophenylacetic acid": "Phg",
    "cyclohexylglycine": "Chg",
    "alpha-cyclohexylglycine": "Chg",
    "cyclopentylglycine": "Gly(cPent)",
    "cyclobutylglycine": "Gly(cBu)",
    "indanylglycine": "Igl",
    "alpha-indanylglycine": "Igl",
    "1-aminocyclopropane-1-carboxylic acid": "ACPrC",
    "1-amino-1-cyclopropanecarboxylic acid": "ACPrC",
    "acc": "ACPrC",
    "1-aminocyclohexane-1-carboxylic acid": "cHex",
    "aminocyclohexane carboxylic acid": "cHex",
    "bicyclo[1.1.1]pentane glycine": "Gly(BCP)",
    "bcp glycine": "Gly(BCP)",
    "bicyclopentylglycine": "Gly(BCP)",
    "difluorocyclohexylalanine": "2F-Cha",
    "4,4-difluorocyclohexylalanine": "2F-Cha",

    # --- 芳香 ---
    "4-fluorophenylalanine": "4-F-Phe",
    "para-fluorophenylalanine": "4-F-Phe",
    "p-fluorophenylalanine": "4-F-Phe",
    "4-chlorophenylalanine": "4-Cl-Phe",
    "para-chlorophenylalanine": "4-Cl-Phe",
    "3-chlorophenylalanine": "3-Cl-Phe",
    "2-fluorophenylalanine": "2-F-Phe",
    "3-chloro-4-fluorophenylalanine": "3-Cl,4-F-Phe",
    "4-cyanophenylalanine": "4-CN-Phe",
    "para-cyanophenylalanine": "4-CN-Phe",
    "4-methylphenylalanine": "4-Me-Phe",
    "3-methylphenylalanine": "3-Me-Phe",
    "4-trifluoromethylphenylalanine": "4-CF3-Phe",
    "para-trifluoromethylphenylalanine": "4-CF3-Phe",
    "biphenylalanine": "Bip",
    "4-biphenylalanine": "Bip",
    "4-phenylphenylalanine": "Bip",
    "3-(4-biphenylyl)alanine": "Bip",
    "2-naphthylalanine": "2-Nal",
    "l-2-naphthylalanine": "2-Nal",
    "beta-naphthylalanine": "2-Nal",
    "3-(2-naphthyl)alanine": "2-Nal",
    "3-(2-naphthyl)-l-alanine": "2-Nal",
    "1-naphthylalanine": "1-Nal",
    "l-1-naphthylalanine": "1-Nal",
    "alpha-naphthylalanine": "1-Nal",
    "3-(1-naphthyl)alanine": "1-Nal",
    "2-pyridylalanine": "2-Pal",
    "3-(2-pyridyl)alanine": "2-Pal",
    "3-pyridylalanine": "3-Pal",
    "3-(3-pyridyl)alanine": "3-Pal",
    "4-pyridylalanine": "4-Pal",
    "3-(4-pyridyl)alanine": "4-Pal",
    "4-methoxyphenylalanine": "4-OMe-Phe",
    "o-methyltyrosine": "4-OMe-Phe",
    "para-methoxyphenylalanine": "4-OMe-Phe",
    "4-aminophenylalanine": "4-NH2-Phe",
    "para-aminophenylalanine": "4-NH2-Phe",
    "3-fluorotyrosine": "3-F-Tyr",
    "3-fluoro-l-tyrosine": "3-F-Tyr",
    "benzothienylalanine": "Bta",
    "beta-benzothienylalanine": "Bta",
    "3-(2-benzothienyl)alanine": "Bta",
    "3-(3-benzothienyl)alanine": "Bta",
    "thienylalanine": "3-Thi",
    "3-thienylalanine": "3-Thi",
    "3-(thiophen-3-yl)alanine": "3-Thi",
    "homophenylalanine": "h-Phe",
    "l-homophenylalanine": "h-Phe",
    "2-amino-4-phenylbutanoic acid": "h-Phe",
    "1,2,3,4-tetrahydroisoquinoline-3-carboxylic acid": "Tic",
    "tetrahydroisoquinolinecarboxylic acid": "Tic",
    "3-methylhistidine": "3-Me-His",
    "tau-methylhistidine": "3-Me-His",
    "4-carbamoylphenylglycine": "4-CONH2-Phg",
    "4-(2-aminoethoxy)phenylalanine": "Phe[4-(2-aminoethoxy)]",
    "3,3-diphenylalanine": "33DPA",
    "alpha,alpha-diphenylalanine": "33DPA",

    # --- Pro 类 ---
    "4-hydroxyproline": "Hyp",
    "trans-4-hydroxyproline": "Hyp",
    "l-4-hydroxyproline": "Hyp",
    "trans-l-hydroxyproline": "Hyp",
    "hydroxyproline": "Hyp",
    "(2s,4r)-4-hydroxyproline": "Hyp",
    "3-hydroxyproline": "Hyp(3-OH)",
    "cis-3-hydroxyproline": "Hyp(3-OH)",
    "4,4-difluoroproline": "Pro(4-F2)",
    "difluoroproline": "Pro(4-F2)",
    "o-ethylhydroxyproline": "Hyp(Et)",
    "o-benzylhydroxyproline": "Hyp(Bzl)",
    "4-oxoproline": "Pro(4-keto)",
    "4-ketoproline": "Pro(4-keto)",
    "4-phenylproline": "Pro(4-Ph)",
    "4-ethylproline": "Pro(4R-Et)",
    "4-propylproline": "Pro(4R-nPr)",
    "homoproline": "hPro",
    "l-homoproline": "hPro",
    "pipecolic acid": "Pip",
    "pipecolinic acid": "Pip",
    "l-pipecolic acid": "Pip",
    "2-piperidinecarboxylic acid": "Pip",
    "homopipecolic acid": "h-Pic",
    "4-oxopipecolic acid": "Pic(4-Oxo)",
    "thiazolidine-4-carboxylic acid": "Thz",
    "thiaproline": "Thz",
    "thioproline": "Thz",
    "l-thiazolidine-4-carboxylic acid": "Thz",
    "azetidine-2-carboxylic acid": "Aze",
    "l-azetidine-2-carboxylic acid": "Aze",
    "3-methylazetidine-2-carboxylic acid": "Aze(3-Me)",
    "3,3-dimethylazetidine-2-carboxylic acid": "Aze(3-Me2)",
    "indoline-2-carboxylic acid": "Idc",
    "l-indoline-2-carboxylic acid": "Idc",
    "octahydroindole-2-carboxylic acid": "Oic",
    "(3as,7as)-octahydroindole-2-carboxylic acid": "Oic",
    "perhydroindole-2-carboxylic acid": "Oic",
    "morpholine-3-carboxylic acid": "Mor",

    # --- Trp 类 ---
    "5-fluorotryptophan": "5-F-Trp",
    "5-fluoro-l-tryptophan": "5-F-Trp",
    "7-fluorotryptophan": "7-F-Trp",
    "7-fluoro-l-tryptophan": "7-F-Trp",
    "5-methyltryptophan": "5-Me-Trp",
    "7-methyltryptophan": "7-Me-Trp",
    "2-methyltryptophan": "2-Me-Trp",
    "methyltryptophan": "Me-Trp",
    "beta-methyltryptophan": "β-Me-Trp",
    "alpha-methyltryptophan": "β-Me-Trp",
    "5-hydroxytryptophan": "5-OH-Trp",
    "5-htp": "5-OH-Trp",
    "oxitriptan": "5-OH-Trp",
    "5-hydroxy-l-tryptophan": "5-OH-Trp",
    "azatryptophan": "AzaTrp",
    "aza-tryptophan": "AzaTrp",
    "dehydrotryptophan": "Dht",
    "2,3-dehydrotryptophan": "Dht",
    "dihydrotryptophan": "Dht",

    # --- Ser/Thr ---
    "o-methylserine": "Ser(Me)",
    "o-methyl-l-serine": "Ser(Me)",
    "o-acetylserine": "Ser(Ac)",
    "o-acetyl-l-serine": "Ser(Ac)",
    "o-benzylserine": "Ser(Bzl)",
    "o-benzyl-l-serine": "Ser(Bzl)",
    "o-isopentylserine": "Ser(iPen)",
    "o-propylserine": "Ser(nPr)",
    "o-n-propylserine": "Ser(nPr)",
    "o-cyclopropylserine": "Ser(cPr)",
    "o-cyclobutylserine": "Ser(cBu)",
    "o-methylthreonine": "Thr(Me)",
    "borono-alanine": "Ala(B(OH)2)",
    "boronoalanine": "Ala(B(OH)2)",
    "3-hydroxy-l-valine": "Val(3-OH)",

    # --- β-AA ---
    "3-aminobutyric acid": "3-Abu",
    "beta-aminobutyric acid": "3-Abu",
    "3-aminobutanoic acid": "3-Abu",
    "3-(trifluoromethyl)-beta-alanine": "3-(CF3)-bAla",
    "3-trifluoromethyl-beta-alanine": "3-(CF3)-bAla",
    "2-aminocyclohexanecarboxylic acid": "2-ACHxC",
    "trans-2-aminocyclohexane carboxylic acid": "2-ACHxC",
    "trans-achc": "2-ACHxC",
    "2-aminocyclopentanecarboxylic acid": "2-ACPnC",
    "trans-2-aminocyclopentane carboxylic acid": "2-ACPnC",
    "alpha-methylphenylalanine": "R-AMPA",
    "r-alpha-methylphenylalanine": "R-AMPA",
    "homohomophenylalanine": "h-Hph",
    "2-amino-5-phenylpentanoic acid": "h-Hph",
}

# ---------------------------------------------------------------------------
# 预编译正则：每个别称 → 单词边界匹配（大小写不敏感）
# 按字符串长度倒序排列，优先匹配最长别称，避免短词吞掉长词
# ---------------------------------------------------------------------------
_SORTED_ALIASES: List[Tuple[re.Pattern, str]] = []

def _build_patterns() -> None:
    seen_canonical: dict = {}
    # 按别称长度降序排列（长优先）
    for alias in sorted(ALIAS_TO_CANONICAL.keys(), key=len, reverse=True):
        canonical = ALIAS_TO_CANONICAL[alias]
        # 转义特殊字符，使用单词边界（支持连字符词如 tert-leucine）
        escaped = re.escape(alias)
        pattern = re.compile(r"(?<![a-zA-Z0-9-])" + escaped + r"(?![a-zA-Z0-9-])",
                             re.IGNORECASE)
        _SORTED_ALIASES.append((pattern, canonical))

_build_patterns()


def tag_compounds(title: str, abstract: str) -> str:
    """
    扫描 title + abstract，返回命中的规范化合物名（逗号分隔）。
    未命中时返回空字符串。
    """
    text = " ".join(filter(None, [str(title or ""), str(abstract or "")]))
    found: list = []
    seen_canonical: set = set()

    for pattern, canonical in _SORTED_ALIASES:
        if canonical in seen_canonical:
            continue
        if pattern.search(text):
            found.append(canonical)
            seen_canonical.add(canonical)

    return ", ".join(found)


def enrich_records(records: list) -> list:
    """就地为每条记录添加 nnaa_compound 字段，返回原列表。"""
    for rec in records:
        rec["nnaa_compound"] = tag_compounds(
            rec.get("title", ""), rec.get("abstract", "")
        )
    return records
