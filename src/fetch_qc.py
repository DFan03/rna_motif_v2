"""Fetch RefSeq/demo RNA records from NCBI and run sequence QC."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import logging
import random
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

LOGGER = logging.getLogger(__name__)

Record = dict[str, Any]

FAMILIES: tuple[str, ...] = (
    "Flaviviridae",
    "Coronaviridae",
    "Picornaviridae",
    "Orthomyxoviridae",
    "Togaviridae",
    "Rhabdoviridae",
    "Caliciviridae",
    "Paramyxoviridae",
    "Retroviridae",
    "Bunyavirales",
    "Nodaviridae",
    "Astroviridae",
)

DEMO_ORGANISMS: dict[str, str] = {
    "Severe acute respiratory syndrome coronavirus 2": "Coronaviridae",
    "Dengue virus 2": "Flaviviridae",
    "Hepatitis C virus": "Flaviviridae",
}

EMAIL = "fandongli1030@gmail.com"
TOOL = "rna_motif_v2"
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

KEPT_COLUMNS: tuple[str, ...] = (
    "accession",
    "set",
    "family",
    "organism",
    "taxid",
    "title",
    "create_date",
    "seq_length",
    "moltype",
    "length",
    "gc",
    "n_masked",
    "alphabet_note",
)

DROP_COLUMNS: tuple[str, ...] = ("accession", "reason")


class EutilsClient:
    """Rate-limited NCBI E-utilities client with retries and timeouts."""

    def __init__(
        self,
        email: str = EMAIL,
        tool: str = TOOL,
        min_interval: float = 0.4,
        timeout: float = 30.0,
        retries: int = 3,
    ) -> None:
        """Initialize client; inputs are NCBI identity, interval seconds, timeout seconds, retries."""
        self.email = email
        self.tool = tool
        self.min_interval = min_interval
        self.timeout = timeout
        self.retries = retries
        self._last_call = 0.0

    def _throttle(self) -> None:
        """Sleep so consecutive HTTP calls are at least min_interval seconds apart."""
        elapsed = time.monotonic() - self._last_call
        if elapsed < self.min_interval:
            time.sleep(self.min_interval - elapsed)

    def get(self, endpoint: str, params: Mapping[str, object]) -> str:
        """Return UTF-8 text from one E-utilities endpoint; inputs endpoint/params, output text."""
        query = {key: str(value) for key, value in params.items()}
        query.update({"email": self.email, "tool": self.tool})
        url = f"{EUTILS_BASE}/{endpoint}?{urllib.parse.urlencode(query)}"

        last_error: Exception | None = None
        for attempt in range(self.retries + 1):
            self._throttle()
            try:
                with urllib.request.urlopen(url, timeout=self.timeout) as response:
                    data = response.read().decode("utf-8")
                return data
            except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, OSError) as exc:
                last_error = exc
                LOGGER.warning(
                    "E-utilities request failed endpoint=%s attempt=%d/%d error=%s",
                    endpoint,
                    attempt + 1,
                    self.retries + 1,
                    exc,
                )
                if attempt < self.retries:
                    time.sleep(min(2 ** attempt, 10))
            finally:
                self._last_call = time.monotonic()

        raise RuntimeError(f"E-utilities request failed after retries: {last_error}")


def _package_versions() -> dict[str, str]:
    """Return installed versions for Python and optional scientific packages."""
    versions = {"python": sys.version.split()[0]}
    for package in ("numpy", "torch", "biopython"):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            continue
    return versions


def _seed_all(seed: int) -> dict[str, str]:
    """Seed random/numpy/torch if available; input integer seed, output package versions."""
    random.seed(seed)
    versions = _package_versions()

    try:
        import numpy as np

        np.random.seed(seed)
    except ImportError:
        pass

    try:
        import torch

        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
    except ImportError:
        pass

    return versions


def _chunks(items: Sequence[str], size: int) -> Iterable[Sequence[str]]:
    """Yield consecutive chunks of size at most size from input sequence."""
    for start in range(0, len(items), size):
        yield items[start : start + size]


def _esearch(client: EutilsClient, term: str, retmax: int) -> list[str]:
    """Return up to retmax nuccore UIDs for one query term."""
    text = client.get(
        "esearch.fcgi",
        {"db": "nuccore", "term": term, "retmax": retmax, "retmode": "json"},
    )
    payload = json.loads(text)
    return [str(uid) for uid in payload.get("esearchresult", {}).get("idlist", [])]


def _esummary(client: EutilsClient, ids: Sequence[str]) -> dict[str, dict[str, Any]]:
    """Return esummary metadata keyed by UID for a batch of nuccore IDs."""
    if not ids:
        return {}
    text = client.get(
        "esummary.fcgi",
        {"db": "nuccore", "id": ",".join(ids), "retmode": "json"},
    )
    payload = json.loads(text)
    result = payload.get("result", {})
    return {
        uid: value
        for uid, value in result.items()
        if uid in ids and isinstance(value, dict)
    }


def _efetch_fasta(client: EutilsClient, ids: Sequence[str]) -> str:
    """Return raw FASTA text from efetch for a batch of nuccore IDs."""
    return client.get(
        "efetch.fcgi",
        {"db": "nuccore", "id": ",".join(ids), "rettype": "fasta", "retmode": "text"},
    )


def _parse_fasta(text: str) -> list[tuple[str, str]]:
    """Parse FASTA text into (header without '>', sequence) pairs."""
    entries: list[tuple[str, str]] = []
    header: str | None = None
    sequence_parts: list[str] = []

    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            if header is not None:
                entries.append((header, "".join(sequence_parts)))
            header = stripped[1:].strip()
            sequence_parts = []
        elif header is not None:
            sequence_parts.append(stripped)

    if header is not None:
        entries.append((header, "".join(sequence_parts)))

    return entries


def _header_accession(header: str) -> str:
    """Return the first whitespace-delimited token from a FASTA header."""
    return header.split(maxsplit=1)[0] if header else ""


def _header_tag(header: str, key: str) -> str:
    """Return value of a bracketed FASTA-header tag like [organism=...], or empty string."""
    marker = f"[{key}="
    start = header.find(marker)
    if start == -1:
        return ""
    start += len(marker)
    end = header.find("]", start)
    if end == -1:
        return ""
    return header[start:end]


def _every_kth(ids: Sequence[str], n: int, seed: int) -> list[str]:
    """Select about n deterministic every-k-th IDs; inputs IDs/n/seed, output selected IDs."""
    if n <= 0 or not ids:
        return []
    if len(ids) <= n:
        return list(ids)

    k = max(1, len(ids) // n)
    rng = random.Random(seed)
    offset = rng.randrange(k) if k > 1 else 0
    chosen = list(ids[offset::k])[:n]

    if len(chosen) < n:
        chosen_set = set(chosen)
        for uid in ids:
            if uid not in chosen_set:
                chosen.append(uid)
                chosen_set.add(uid)
                if len(chosen) >= n:
                    break

    return chosen


def _make_record(
    accession: str,
    header: str,
    sequence: str,
    summary: Mapping[str, Any],
    family: str,
    set_name: str,
    organism_hint: str | None = None,
) -> Record:
    """Build one metadata record from summary/FASTA fields; output record dict."""
    title = str(summary.get("title") or "")
    if not title:
        parts = header.split(maxsplit=1)
        title = parts[1] if len(parts) > 1 else ""

    organism = str(summary.get("organism") or organism_hint or _header_tag(header, "organism"))
    taxid = summary.get("taxid")
    if taxid in (None, ""):
        taxid = _header_tag(header, "taxid")

    seq_length_raw = summary.get("slen")
    try:
        seq_length = int(seq_length_raw) if seq_length_raw not in (None, "") else len(sequence)
    except (TypeError, ValueError):
        seq_length = len(sequence)

    moltype = str(summary.get("moltype") or summary.get("biomol") or "")
    create_date = str(summary.get("createdate") or summary.get("create_date") or "")

    return {
        "accession": accession,
        "set": set_name,
        "family": family,
        "organism": organism,
        "taxid": taxid,
        "title": title,
        "create_date": create_date,
        "seq_length": seq_length,
        "moltype": moltype,
        "sequence": sequence,
        "alphabet_note": "DNA alphabet (ACGT); T retained (GenBank cDNA convention)",
    }


def _fetch_ids(
    client: EutilsClient,
    ids: Sequence[str],
    family: str,
    set_name: str,
    organism_hint: str | None = None,
) -> tuple[list[Record], list[str]]:
    """Fetch summaries and FASTA for IDs in batches of 50; output records and raw text chunks."""
    records: list[Record] = []
    raw_texts: list[str] = []
    if not ids:
        return records, raw_texts

    batches = list(_chunks(ids, 50))
    start_time = time.monotonic()

    for batch_index, batch in enumerate(batches, start=1):
        summaries = _esummary(client, batch)
        fasta_text = _efetch_fasta(client, batch)
        raw_texts.append(fasta_text)
        entries = _parse_fasta(fasta_text)
        summary_by_accession = {
            str(summary.get("accessionversion", "")): summary
            for summary in summaries.values()
        }

        for entry_index, (header, sequence) in enumerate(entries):
            accession = _header_accession(header)
            summary: Mapping[str, Any] = {}
            if entry_index < len(batch):
                summary = summaries.get(batch[entry_index], {})

            summary_accession = str(summary.get("accessionversion", "")) if summary else ""
            if not summary or (summary_accession and summary_accession != accession):
                summary = summary_by_accession.get(accession, summary)

            records.append(
                _make_record(
                    accession=accession,
                    header=header,
                    sequence=sequence,
                    summary=summary,
                    family=family,
                    set_name=set_name,
                    organism_hint=organism_hint,
                )
            )

        elapsed = time.monotonic() - start_time
        percent = 100.0 * batch_index / len(batches)
        LOGGER.info(
            "fetched set=%s family=%s batch=%d/%d done=%d/%d percent=%.1f elapsed=%.1fs",
            set_name,
            family,
            batch_index,
            len(batches),
            batch_index * 50,
            len(ids),
            percent,
            elapsed,
        )

    return records, raw_texts


def _fetch_refseq(
    client: EutilsClient,
    families: Sequence[str],
    per_family: int,
) -> tuple[list[Record], list[str], list[dict[str, Any]]]:
    """Fetch first per_family RefSeq RNA records per family; output records, raw text, queries."""
    records: list[Record] = []
    raw_texts: list[str] = []
    queries: list[dict[str, Any]] = []

    for family in families:
        term = (
            f'"{family}"[Organism] AND srcdb_refseq[PROP] '
            "AND 1000:35000[SLEN]"
        )
        queries.append(
            {"set": "refseq", "family": family, "retmax": per_family, "term": term}
        )
        ids = _esearch(client, term, per_family)
        LOGGER.info("esearch refseq family=%s found=%d requested=%d", family, len(ids), per_family)
        family_records, family_raw = _fetch_ids(client, ids, family, "refseq")
        records.extend(family_records)
        raw_texts.extend(family_raw)

    return records, raw_texts, queries


def _fetch_demo(
    client: EutilsClient,
    demo_n: int,
    seed: int,
) -> tuple[list[Record], list[str], list[dict[str, Any]]]:
    """Fetch deterministic every-k-th demo GenBank records; output records, raw text, queries."""
    records: list[Record] = []
    raw_texts: list[str] = []
    queries: list[dict[str, Any]] = []

    for organism, family in DEMO_ORGANISMS.items():
        term = (
            f'"{organism}"[Organism] AND biomol_genomic[PROP] '
            "AND 7000:35000[SLEN] NOT srcdb_refseq[PROP]"
        )
        queries.append(
            {"set": "demo", "family": family, "retmax": 500, "term": term}
        )
        ids = _esearch(client, term, 500)
        selected = _every_kth(ids, demo_n, seed)
        LOGGER.info(
            "esearch demo organism=%s found=%d selected=%d requested=%d",
            organism,
            len(ids),
            len(selected),
            demo_n,
        )
        organism_records, organism_raw = _fetch_ids(
            client,
            selected,
            family,
            "demo",
            organism_hint=organism,
        )
        records.extend(organism_records)
        raw_texts.extend(organism_raw)

    return records, raw_texts, queries


def qc_records(
    records: Sequence[Mapping[str, Any]],
    min_len: int = 200,
    max_len: int = 35000,
    max_nonacgt_frac: float = 0.01,
) -> tuple[list[Record], list[dict[str, str]]]:
    """QC records; inputs records and length/non-ACGT thresholds, outputs (kept, drops)."""
    kept: list[Record] = []
    drops: list[dict[str, str]] = []
    seen_sequences: set[str] = set()

    for original in records:
        record = dict(original)
        accession = str(record.get("accession", ""))
        sequence = "".join(str(record.get("sequence", "") or "").split()).upper()

        if not sequence:
            drops.append({"accession": accession, "reason": "empty_sequence"})
            continue

        length = len(sequence)
        if length < min_len:
            drops.append({"accession": accession, "reason": "too_short"})
            continue
        if length > max_len:
            drops.append({"accession": accession, "reason": "too_long"})
            continue

        taxid = record.get("taxid")
        taxid_text = "" if taxid is None else str(taxid).strip()
        organism = str(record.get("organism", "") or "").strip()
        if not organism or taxid_text in {"", "0", "None", "nan"}:
            drops.append({"accession": accession, "reason": "missing_taxonomy"})
            continue

        nonstandard_count = sum(base not in {"A", "C", "G", "T"} for base in sequence)
        nonstandard_frac = nonstandard_count / length
        if nonstandard_frac > max_nonacgt_frac:
            drops.append({"accession": accession, "reason": "nonstandard_frac>1%"})
            continue

        digest = hashlib.sha1(sequence.encode("utf-8")).hexdigest()
        if digest in seen_sequences:
            drops.append({"accession": accession, "reason": "exact_duplicate_sequence"})
            continue
        seen_sequences.add(digest)

        masked_sequence = "".join(base if base in {"A", "C", "G", "T"} else "N" for base in sequence)
        gc = (masked_sequence.count("G") + masked_sequence.count("C")) / length

        record["sequence"] = masked_sequence
        record["length"] = length
        record["gc"] = gc
        record["n_masked"] = nonstandard_count
        kept.append(record)

    return kept, drops


def _write_raw_fasta(path: Path, raw_texts: Sequence[str]) -> None:
    """Write concatenated raw efetch FASTA text to path; no sequence modification."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for text in raw_texts:
            handle.write(text)
            if text and not text.endswith("\n"):
                handle.write("\n")


def _write_fasta(path: Path, records: Sequence[Mapping[str, Any]]) -> None:
    """Write records as wrapped FASTA to path using accession/family/set headers."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for record in records:
            accession = str(record.get("accession", ""))
            family = str(record.get("family", ""))
            set_name = str(record.get("set", ""))
            sequence = str(record.get("sequence", ""))
            handle.write(f">{accession} family={family} set={set_name}\n")
            for start in range(0, len(sequence), 70):
                handle.write(sequence[start : start + 70] + "\n")


def _write_tsv(
    path: Path,
    rows: Sequence[Mapping[str, Any]],
    columns: Sequence[str],
) -> None:
    """Write mapping rows to a tab-separated file with the given ordered columns."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(columns),
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def _find_key(payload: Any, key: str) -> Any:
    """Recursively search nested JSON-like payload for key; output first match or None."""
    if isinstance(payload, dict):
        for current_key, value in payload.items():
            if str(current_key).lower() == key.lower():
                return value
            found = _find_key(value, key)
            if found is not None:
                return found
    elif isinstance(payload, list):
        for value in payload:
            found = _find_key(value, key)
            if found is not None:
                return found
    return None


def _fetch_dbbuild(client: EutilsClient) -> str:
    """Return NCBI nuccore dbbuild string from einfo, or 'unknown' on failure."""
    try:
        text = client.get("einfo.fcgi", {"db": "nuccore", "retmode": "json"})
        payload = json.loads(text)
        dbbuild = _find_key(payload, "dbbuild")
        return str(dbbuild) if dbbuild else "unknown"
    except (RuntimeError, json.JSONDecodeError, ValueError) as exc:
        LOGGER.warning("could not fetch nuccore dbbuild: %s", exc)
        return "unknown"


def _write_provenance(
    path: Path,
    queries: Sequence[Mapping[str, Any]],
    dbbuild: str,
    fetch_date: str,
    args: argparse.Namespace,
    versions: Mapping[str, str],
) -> None:
    """Write PROVENANCE.md with fetch date, dbbuild, parameters, versions, and exact queries."""
    lines: list[str] = [
        "# Provenance",
        "",
        f"- Fetch date (UTC): {fetch_date}",
        f"- NCBI nuccore dbbuild: {dbbuild}",
        f"- Seed: {args.seed}",
        f"- Output directory: {args.out}",
        f"- Dry run: {args.dry}",
        "",
        "## Parameters",
        "",
        f"- per_family: {args.per_family}",
        f"- demo_n: {args.demo_n}",
        f"- min_len: 200",
        f"- max_len: 35000",
        f"- max_nonacgt_frac: 0.01",
        "",
        "## Package versions",
        "",
    ]
    for package, version in sorted(versions.items()):
        lines.append(f"- {package}: {version}")

    lines.extend(["", "## Queries", ""])
    for query in queries:
        lines.append(
            f"- set={query['set']} family={query['family']} "
            f"retmax={query['retmax']} term: `{query['term']}`"
        )

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _log_summary(
    kept: Sequence[Mapping[str, Any]],
    drops: Sequence[Mapping[str, Any]],
    all_records: Sequence[Mapping[str, Any]],
) -> None:
    """Log kept/drop counts per family and per drop reason."""
    family_totals = Counter(str(record.get("family", "")) for record in all_records)
    family_kept = Counter(str(record.get("family", "")) for record in kept)

    for family in sorted(family_totals):
        LOGGER.info(
            "summary family=%s kept=%d total=%d",
            family,
            family_kept.get(family, 0),
            family_totals[family],
        )

    reasons = Counter(str(drop.get("reason", "")) for drop in drops)
    for reason in sorted(reasons):
        LOGGER.info("summary drop_reason=%s count=%d", reason, reasons[reason])


def main(argv: Sequence[str] | None = None) -> int:
    """Run CLI fetch/QC pipeline; inputs argv, output process exit code."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("data"), help="Output data directory")
    parser.add_argument("--per-family", type=int, default=40, help="RefSeq records per family")
    parser.add_argument("--demo-n", type=int, default=50, help="Demo records per organism")
    parser.add_argument(
        "--dry",
        type=int,
        default=None,
        help="Tiny run: 2 families, N records each, no demo set",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        stream=sys.stderr,
    )

    if args.per_family < 1:
        parser.error("--per-family must be >= 1")
    if args.demo_n < 1 and args.dry is None:
        parser.error("--demo-n must be >= 1")
    if args.dry is not None and args.dry < 1:
        parser.error("--dry must be >= 1")

    versions = _seed_all(args.seed)
    out_dir = Path(args.out)
    raw_dir = out_dir / "raw"
    qc_dir = out_dir / "qc"
    raw_dir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)

    LOGGER.info(
        "start out=%s per_family=%d demo_n=%d dry=%s seed=%d",
        out_dir,
        args.per_family,
        args.demo_n,
        args.dry,
        args.seed,
    )

    if args.dry is not None:
        families = FAMILIES[:2]
        per_family = args.dry
        include_demo = False
    else:
        families = FAMILIES
        per_family = args.per_family
        include_demo = True

    fetch_date = datetime.now(timezone.utc).isoformat(timespec="seconds")
    client = EutilsClient()
    dbbuild = _fetch_dbbuild(client)
    queries: list[dict[str, Any]] = []

    refseq_records, refseq_raw, refseq_queries = _fetch_refseq(client, families, per_family)
    queries.extend(refseq_queries)
    _write_raw_fasta(raw_dir / "refseq.fa", refseq_raw)

    demo_records: list[Record] = []
    demo_raw: list[str] = []
    if include_demo:
        demo_records, demo_raw, demo_queries = _fetch_demo(client, args.demo_n, args.seed)
        queries.extend(demo_queries)
        _write_raw_fasta(raw_dir / "demo.fa", demo_raw)

    refseq_kept, refseq_drops = qc_records(refseq_records)
    if include_demo:
        demo_kept, demo_drops = qc_records(demo_records)
    else:
        demo_kept, demo_drops = [], []

    kept = refseq_kept + demo_kept
    drops = refseq_drops + demo_drops
    all_records = refseq_records + demo_records

    _write_fasta(qc_dir / "kept_refseq.fa", refseq_kept)
    if include_demo:
        _write_fasta(qc_dir / "kept_demo.fa", demo_kept)

    _write_tsv(qc_dir / "records.tsv", kept, KEPT_COLUMNS)
    _write_tsv(qc_dir / "drops.tsv", drops, DROP_COLUMNS)
    _write_provenance(out_dir / "PROVENANCE.md", queries, dbbuild, fetch_date, args, versions)

    _log_summary(kept, drops, all_records)
    LOGGER.info("done output_dir=%s kept=%d dropped=%d", out_dir, len(kept), len(drops))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
