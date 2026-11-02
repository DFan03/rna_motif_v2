# Provenance

- Fetch date (UTC): 2026-11-02T02:02:00+00:00
- NCBI nuccore dbbuild: Build261031-0025m.1
- Seed: 42
- Output directory: data
- Dry run: None

## Parameters

- per_family: 40
- demo_n: 50
- min_len: 200
- max_len: 35000
- max_nonacgt_frac: 0.01

## Package versions

- biopython: 1.84
- numpy: 1.26.4
- python: 3.12.7

## Queries

- set=refseq family=Flaviviridae retmax=40 term: `"Flaviviridae"[Organism] AND srcdb_refseq[PROP] AND 1000:35000[SLEN]`
- set=refseq family=Coronaviridae retmax=40 term: `"Coronaviridae"[Organism] AND srcdb_refseq[PROP] AND 1000:35000[SLEN]`
- set=refseq family=Picornaviridae retmax=40 term: `"Picornaviridae"[Organism] AND srcdb_refseq[PROP] AND 1000:35000[SLEN]`
- set=refseq family=Orthomyxoviridae retmax=40 term: `"Orthomyxoviridae"[Organism] AND srcdb_refseq[PROP] AND 1000:35000[SLEN]`
- set=refseq family=Togaviridae retmax=40 term: `"Togaviridae"[Organism] AND srcdb_refseq[PROP] AND 1000:35000[SLEN]`
- set=refseq family=Rhabdoviridae retmax=40 term: `"Rhabdoviridae"[Organism] AND srcdb_refseq[PROP] AND 1000:35000[SLEN]`
- set=refseq family=Caliciviridae retmax=40 term: `"Caliciviridae"[Organism] AND srcdb_refseq[PROP] AND 1000:35000[SLEN]`
- set=refseq family=Paramyxoviridae retmax=40 term: `"Paramyxoviridae"[Organism] AND srcdb_refseq[PROP] AND 1000:35000[SLEN]`
- set=refseq family=Retroviridae retmax=40 term: `"Retroviridae"[Organism] AND srcdb_refseq[PROP] AND 1000:35000[SLEN]`
- set=refseq family=Bunyavirales retmax=40 term: `"Bunyavirales"[Organism] AND srcdb_refseq[PROP] AND 1000:35000[SLEN]`
- set=refseq family=Nodaviridae retmax=40 term: `"Nodaviridae"[Organism] AND srcdb_refseq[PROP] AND 1000:35000[SLEN]`
- set=refseq family=Astroviridae retmax=40 term: `"Astroviridae"[Organism] AND srcdb_refseq[PROP] AND 1000:35000[SLEN]`
- set=demo family=Coronaviridae retmax=500 term: `"Severe acute respiratory syndrome coronavirus 2"[Organism] AND biomol_genomic[PROP] AND 7000:35000[SLEN] NOT srcdb_refseq[PROP]`
- set=demo family=Flaviviridae retmax=500 term: `"Dengue virus 2"[Organism] AND biomol_genomic[PROP] AND 7000:35000[SLEN] NOT srcdb_refseq[PROP]`
- set=demo family=Flaviviridae retmax=500 term: `"Hepatitis C virus"[Organism] AND biomol_genomic[PROP] AND 7000:35000[SLEN] NOT srcdb_refseq[PROP]`
