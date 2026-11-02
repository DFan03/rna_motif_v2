"""Unit tests for qc_records only; no network calls."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fetch_qc import qc_records


Record = dict[str, Any]


@pytest.fixture()
def base_record() -> Record:
    """Return one tiny synthetic normal record for QC tests."""
    return {
        "accession": "TEST.1",
        "set": "refseq",
        "family": "Flaviviridae",
        "organism": "Test virus",
        "taxid": "12345",
        "title": "synthetic test record",
        "create_date": "2024-01-01",
        "seq_length": 300,
        "moltype": "rna",
        "sequence": "ACGT" * 75,
        "alphabet_note": "DNA alphabet (ACGT); T retained (GenBank cDNA convention)",
    }


def _record(base: Record, **updates: Any) -> Record:
    """Return a copy of base_record with selected fields replaced."""
    record = dict(base)
    record.update(updates)
    return record


def test_normal_record_kept(base_record: Record) -> None:
    """Normal 300 bp ACGT record is kept with expected GC and no masking."""
    kept, drops = qc_records([base_record])

    assert not drops
    assert len(kept) == 1
    assert kept[0]["sequence"] == "ACGT" * 75
    assert kept[0]["length"] == 300
    assert kept[0]["gc"] == pytest.approx(0.5)
    assert kept[0]["n_masked"] == 0


def test_empty_records() -> None:
    """Empty input returns empty kept and drops lists."""
    kept, drops = qc_records([])
    assert kept == []
    assert drops == []


def test_short_record_dropped(base_record: Record) -> None:
    """Record shorter than min_len is dropped as too_short."""
    record = _record(base_record, sequence="ACGT" * 2)
    kept, drops = qc_records([record])

    assert not kept
    assert drops == [{"accession": "TEST.1", "reason": "too_short"}]


def test_long_record_dropped(base_record: Record) -> None:
    """Record longer than max_len is dropped as too_long."""
    record = _record(base_record, sequence="A" * 35001)
    kept, drops = qc_records([record])

    assert not kept
    assert drops == [{"accession": "TEST.1", "reason": "too_long"}]


def test_empty_sequence_dropped(base_record: Record) -> None:
    """Record with empty sequence is dropped as empty_sequence."""
    record = _record(base_record, sequence="")
    kept, drops = qc_records([record])

    assert not kept
    assert drops == [{"accession": "TEST.1", "reason": "empty_sequence"}]


def test_nonstandard_fraction_dropped(base_record: Record) -> None:
    """Record with more than 1% non-ACGT bases is dropped."""
    record = _record(base_record, sequence="N" * 5 + "A" * 295)
    kept, drops = qc_records([record])

    assert not kept
    assert drops == [{"accession": "TEST.1", "reason": "nonstandard_frac>1%"}]


def test_iupac_masking_at_threshold(base_record: Record) -> None:
    """Record with exactly 1% IUPAC/non-ACGT bases is kept and masked to N."""
    sequence_list = ["A"] * 300
    sequence_list[0] = "N"
    sequence_list[1] = "R"
    sequence_list[2] = "Y"
    record = _record(base_record, sequence="".join(sequence_list))

    kept, drops = qc_records([record])

    assert not drops
    assert len(kept) == 1
    assert kept[0]["n_masked"] == 3
    assert kept[0]["sequence"][:3] == "NNN"
    assert kept[0]["sequence"][3:] == "A" * 297


def test_exact_duplicate_sequence_dropped(base_record: Record) -> None:
    """Second record with identical uppercase sequence is dropped as duplicate."""
    first = _record(base_record, accession="FIRST.1")
    second = _record(base_record, accession="SECOND.1", sequence=("ACGT" * 75).lower())

    kept, drops = qc_records([first, second])

    assert len(kept) == 1
    assert kept[0]["accession"] == "FIRST.1"
    assert drops == [{"accession": "SECOND.1", "reason": "exact_duplicate_sequence"}]


def test_missing_taxonomy_dropped(base_record: Record) -> None:
    """Record missing taxid is dropped as missing_taxonomy."""
    record = _record(base_record, taxid=None)
    kept, drops = qc_records([record])

    assert not kept
    assert drops == [{"accession": "TEST.1", "reason": "missing_taxonomy"}]


def test_missing_organism_dropped(base_record: Record) -> None:
    """Record missing organism is dropped as missing_taxonomy."""
    record = _record(base_record, organism="")
    kept, drops = qc_records([record])

    assert not kept
    assert drops == [{"accession": "TEST.1", "reason": "missing_taxonomy"}]
