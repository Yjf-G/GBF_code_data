# Zenodo deposit guide

This file documents how the Zenodo record for this repository is created and
which metadata is used. It is the operational companion of `CITATION.cff` and
`.zenodo.json` (the authoritative machine-readable metadata). Zenodo's manual
web upload does **not** read those two files automatically — paste the values
below into the web form (the GitHub–Zenodo integration would read them
automatically; see "Alternative route" at the end).

Deposit title: **Code and data for: Urgent and balanced action needed to
achieve global biodiversity framework** — one record, Zenodo type *Software*,
license *MIT*, version *1.0.0*.

---

## 1. The archive

The ready-made deposit archive is **`GBF_code_data_1.0.0.zip`** (in this
folder; ~63 MB). It contains every file of this repository under a single
top-level folder `GBF_code_data_1.0.0/`, so the archive is self-describing
(`LICENSE`, `CITATION.cff`, `.zenodo.json` travel with it) and does not
include itself. The zip is excluded from git (`.gitignore`), so the GitHub
repository does not carry a duplicate.

To rebuild the archive after any change, run this from the folder that
contains `upload/` (or adjust the paths):

```sh
python -I -c "
import os, zipfile
SRC='upload'; OUT='GBF_code_data_1.0.0.zip'; WRAP='GBF_code_data_1.0.0/'
with zipfile.ZipFile(OUT,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for dp,dn,fns in os.walk(SRC):
        dn.sort()
        for fn in sorted(fns):
            if fn==OUT or fn.startswith('upload.7z'): continue
            p=os.path.join(dp,fn)
            z.write(p, WRAP+os.path.relpath(p,SRC).replace(os.sep,'/'))
"
mv GBF_code_data_1.0.0.zip upload/
```

The repository is ~150 MB uncompressed — far below the Zenodo limits
(50 GB per file).

## 2. Fill the Zenodo web form

On <https://zenodo.org> → **New upload** → upload `GBF_code_data_1.0.0.zip` →
fill in:

| Form field | Value |
|---|---|
| Resource type | **Software** |
| Title | `Code and data for: Urgent and balanced action needed to achieve global biodiversity framework` |
| Creators | In this order (copy names and affiliations from `.zenodo.json`, one entry per author): Yi, Jianfeng; Gregory, Richard D.; Cao, Mingchang; Cao, Yun; Dubois, Grégoire; Gonzalez, Andrew; Lin, Hongmei; Liu, Wei; Delli, Giacomo; Wu, Jianyong; Zhang, Riquan; Yu, Dandan; Xu, Haigen |
| Description | Paste the three-paragraph description from `.zenodo.json` (plain text is fine) |
| License | **MIT** (MIT License) |
| Version | `1.0.0` |
| Keywords | Kunming–Montreal Global Biodiversity Framework; GBF; biodiversity; conservation targets; indicators; achievement rate; projection; reproducibility |
| Related works | *Is supplement to* → the manuscript DOI (add once available; can be edited in later without changing this record's DOI) |
| Funding / Communities / References | leave empty unless known |

Then **Save** (the record stays a draft with a reserved DOI).

## 3. Review and publish

- A **private link** can be generated for the draft (right-hand panel) and sent
  to editors/reviewers before publication — the reserved DOI does not resolve
  publicly until the record is published.
- Publish before manuscript submission so that the DOI cited in the Code
  availability statement and in the reference list is live.

## 4. After publication

1. Put the Zenodo DOI into `CITATION.cff` (`doi:` field) and, optionally, the
   `preferred-citation` block once the manuscript has a DOI.
2. Add the citation to the "Citation" section of `README.md`.
3. Cite the record in the manuscript:
   - *Code availability* — "All analysis code used in this study is publicly
     available in Zenodo at https://doi.org/10.5281/zenodo.XXXXXXX (ref. N)."
   - Reference list (Nature style) —
     `Yi, J. et al. Code and data for: Urgent and balanced action needed to achieve global biodiversity framework. Zenodo https://doi.org/10.5281/zenodo.XXXXXXX (2026).`
   - If the derived data in `data/` are cited in *Data availability*, the same
     record and DOI can be referenced there.
4. After acceptance, edit the record metadata and add the manuscript DOI under
   **Related works** (*Is supplement to*) — metadata edits keep the same DOI.

## Alternative route: GitHub–Zenodo integration

If this repository is pushed to GitHub, creating a release there and connecting
the repository on <https://zenodo.org> ("GitHub" in the account menu) archives
the release automatically, mints the DOI, and reads `.zenodo.json` /
`CITATION.cff` for the metadata. The reserved DOI is written back to the
repository automatically.
