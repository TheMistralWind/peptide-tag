"""Find a peptide sequence inside the real human proteome.

This is the honest replacement for the invented "therapeutic applications"
section. Every hit is a literal substring match against reviewed human
Swiss-Prot entries, so anything we say here is checkable.

UniProt's own peptide-search endpoints (peptidesearch.uniprot.org and the
PIR peptidematch API) were timing out and returning 503 when this was built,
so we ship the proteome and search it in-process instead. It is 7.2 MB
gzipped, loads in well under a second, and answers in single-digit
milliseconds with no network call and nothing that can go down.

Refresh with: python -m proteome --refresh
"""

from __future__ import annotations

import bisect
import gzip
import os
import re
import threading
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

DATA = Path(__file__).resolve().parent / "data" / "human_swissprot.fasta.gz"
SOURCE_URL = (
    "https://rest.uniprot.org/uniprotkb/stream"
    "?query=(organism_id:9606)+AND+(reviewed:true)&format=fasta&compressed=true"
)

# Sequences are joined with a character that cannot appear in a residue run,
# so no match can straddle two proteins.
_SEPARATOR = "|"

_HEADER = re.compile(
    r"^>\w\w\|(?P<acc>[^|]+)\|(?P<entry>\S+)\s+(?P<desc>.*?)(?:\s+OS=|$)"
)
_GENE = re.compile(r"\bGN=(\S+)")


@dataclass(frozen=True)
class Protein:
    accession: str
    entry: str
    name: str
    gene: str
    length: int


@dataclass(frozen=True)
class Hit:
    protein: Protein
    position: int  # 1-based residue position of the match

    @property
    def url(self) -> str:
        return f"https://www.uniprot.org/uniprotkb/{self.protein.accession}"

    @property
    def blurb(self) -> str:
        """Why a non-biologist should care, when we can say."""
        return FAMOUS.get(self.protein.gene, "")

    @property
    def famous(self) -> bool:
        return self.protein.gene in FAMOUS


# Short names hit hundreds of proteins, and the first match in file order is
# almost always something like "Probable oxidoreductase PXDNL". Ranking by
# recognisability is the difference between a fact and a delight. The blurb is
# what makes the hit land, so only proteins worth a sentence are in here.
FAMOUS: dict[str, str] = {
    "TTN": "the largest protein in your body",
    "TP53": "the tumour suppressor, the most studied gene in medicine",
    "INS": "insulin",
    "ALB": "albumin, the most abundant protein in your blood",
    "HBB": "haemoglobin, which carries oxygen around your body",
    "ACTB": "actin, part of the scaffolding inside every one of your cells",
    "COL1A1": "collagen, which is most of what holds you together",
    "KRT1": "keratin, the stuff of hair and fingernails",
    "DMD": "dystrophin, the protein missing in muscular dystrophy",
    "APP": "the amyloid precursor protein, central to Alzheimer's research",
    "SNCA": "alpha-synuclein, the protein that clumps in Parkinson's",
    "MAPT": "tau, the other protein implicated in dementia",
    "PRNP": "the prion protein, the one behind mad cow disease",
    "CFTR": "the channel that fails in cystic fibrosis",
    "MB": "myoglobin, which stores oxygen in your muscles",
    "SOD1": "superoxide dismutase, one of your antioxidant defences",
    "OXT": "oxytocin, the bonding hormone",
    "POMC": "the precursor your body cuts up to make endorphins",
    "GCG": "glucagon, insulin's opposite number",
    "LEP": "leptin, the hormone that tells you that you are full",
    "EGFR": "the growth factor receptor targeted by many cancer drugs",
    "BRCA1": "the breast cancer susceptibility protein",
    "BRCA2": "the second breast cancer susceptibility protein",
    "APOE": "the cholesterol carrier whose variants affect Alzheimer's risk",
    "VWF": "von Willebrand factor, which helps your blood clot",
    "FGB": "fibrinogen, which becomes the mesh of a scab",
    "GAPDH": "the workhorse enzyme of sugar metabolism",
    "MYC": "one of the most powerful cancer-driving genes known",
    "KRAS": "the switch mutated in a third of all cancers",
    "PTEN": "a tumour suppressor that brakes cell growth",
    "MTOR": "the master regulator of cell growth and ageing",
    "ACE2": "the receptor SARS-CoV-2 uses to get into your cells",
    "CRP": "C-reactive protein, what a blood test measures for inflammation",
    "TG": "thyroglobulin, the raw material for thyroid hormone",
    "IGF1": "insulin-like growth factor, which made you grow",
    "GH1": "growth hormone",
    "VIM": "vimentin, part of your cells' internal skeleton",
    "LMNA": "lamin, the lining of the nucleus, tied to premature ageing",
    "CAT": "catalase, which destroys hydrogen peroxide inside you",
    "F8": "clotting factor VIII, the one missing in haemophilia",
    "MYH7": "the motor protein that contracts your heart",
    "FN1": "fibronectin, the glue between your cells",
    "SERPINA1": "alpha-1 antitrypsin, which protects your lungs",
    "HTT": "huntingtin, the protein behind Huntington's disease",
}

_VAGUE = ("uncharacterized", "probable", "putative", "unknown")


def _notability(protein: "Protein") -> tuple:
    """Sort key: famous first, then named, then anything. Higher is better.

    Length is negated so that among equally recognisable proteins the smaller
    one wins. Landing inside a 100-residue hormone is a better story than
    landing inside a 9,000-residue giant, where almost any short run matches.
    """
    if protein.gene in FAMOUS:
        return (2, -protein.length)
    lowered = protein.name.lower()
    if any(w in lowered for w in _VAGUE):
        return (0, -protein.length)
    return (1, -protein.length)


@dataclass(frozen=True)
class Fragment:
    text: str
    start: int  # 0-based offset into the name's sequence
    whole: bool  # True when the entire name matched
    hits: list["Hit"]
    occurrences: int

    @property
    def end(self) -> int:
        return self.start + len(self.text)


class Proteome:
    def __init__(self, path: Path = DATA):
        self.path = path
        self.blob = ""
        self.proteins: list[Protein] = []
        self.starts: list[int] = []
        self.residues = 0
        self.loaded = False

    def load(self) -> None:
        if self.loaded or not self.path.exists():
            return
        chunks: list[str] = []
        proteins: list[Protein] = []
        starts: list[int] = []
        offset = 0
        header: str | None = None
        parts: list[str] = []

        def flush() -> None:
            nonlocal offset, header, parts
            if header is None:
                return
            seq = "".join(parts)
            m = _HEADER.match(header)
            if m and seq:
                gene = _GENE.search(header)
                proteins.append(Protein(
                    accession=m.group("acc"),
                    entry=m.group("entry"),
                    name=m.group("desc").strip(),
                    gene=gene.group(1) if gene else "",
                    length=len(seq),
                ))
                starts.append(offset)
                chunks.append(seq)
                chunks.append(_SEPARATOR)
                offset += len(seq) + len(_SEPARATOR)
            header, parts = None, []

        with gzip.open(self.path, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    flush()
                    header = line.rstrip("\n")
                else:
                    parts.append(line.strip())
            flush()

        self.blob = "".join(chunks)
        self.proteins = proteins
        self.starts = starts
        self.residues = sum(p.length for p in proteins)
        self.loaded = True

    def find(self, sequence: str, limit: int = 5, scan: int = 500) -> list[Hit]:
        """Human proteins containing this sequence, most recognisable first.

        A four-letter name matches hundreds of proteins, so returning them in
        file order buries titin under a list of uncharacterised ORFs. We scan up
        to `scan` occurrences and rank what we find.
        """
        if not self.loaded or len(sequence) < 3 or not sequence:
            return []
        found: list[Hit] = []
        seen: set[str] = set()
        start = 0
        while len(found) < scan:
            idx = self.blob.find(sequence, start)
            if idx == -1:
                break
            start = idx + 1
            i = bisect.bisect_right(self.starts, idx) - 1
            if i < 0:
                continue
            protein = self.proteins[i]
            if protein.accession in seen:
                continue
            seen.add(protein.accession)
            found.append(Hit(protein=protein, position=idx - self.starts[i] + 1))
        found.sort(key=lambda h: _notability(h.protein), reverse=True)
        return found[:limit]

    def count(self, sequence: str) -> int:
        """Total occurrences, used to tell 'rare' from 'everywhere'."""
        if not self.loaded or len(sequence) < 3:
            return 0
        return self.blob.count(sequence)

    def longest_fragment(self, sequence: str, min_len: int = 3) -> Fragment | None:
        """The longest run of the name that exists inside a real human protein.

        A whole-name match only happens for about a third of names, and
        essentially never above six letters, so asking "how much of your name is
        in there" gives everybody a result and makes longer names interesting
        instead of hopeless. Ties go to the fragment nearest the start of the
        name, which reads better.
        """
        if not self.loaded or len(sequence) < min_len:
            return None
        for length in range(len(sequence), min_len - 1, -1):
            for start in range(0, len(sequence) - length + 1):
                fragment = sequence[start:start + length]
                if self.blob.find(fragment) != -1:
                    return Fragment(
                        text=fragment,
                        start=start,
                        whole=(length == len(sequence)),
                        hits=self.find(fragment, limit=4),
                        occurrences=self.count(fragment),
                    )
        return None


_proteome = Proteome()
_lock = threading.Lock()


def get() -> Proteome:
    if not _proteome.loaded:
        with _lock:
            if not _proteome.loaded:
                _proteome.load()
    return _proteome


def search(sequence: str, limit: int = 5) -> list[Hit]:
    return get().find(sequence, limit=limit)


@lru_cache(maxsize=4096)
def longest_fragment(sequence: str) -> Fragment | None:
    """Cached: names repeat, and a miss costs a scan per candidate length."""
    return get().longest_fragment(sequence)


def _refresh() -> None:
    import urllib.request
    DATA.parent.mkdir(parents=True, exist_ok=True)
    tmp = DATA.with_suffix(".tmp")
    with urllib.request.urlopen(SOURCE_URL) as r, open(tmp, "wb") as out:
        out.write(r.read())
    os.replace(tmp, DATA)
    print(f"wrote {DATA} ({DATA.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    import sys
    if "--refresh" in sys.argv:
        _refresh()
    p = get()
    print(f"{len(p.proteins):,} proteins, {p.residues:,} residues")
