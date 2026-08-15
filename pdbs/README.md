## PDBs

Reference structures used for local docking (not committed — `*.pdb` files are gitignored).

- `3PTB.pdb` — full trypsin crystal structure, used by `calibrate.py` (positive-control
  benzamidine docking check) and as the rigid-crystal baseline receptor in
  `ensemble_auditor.py`. Download it with:
  ```bash
  curl -L https://files.rcsb.org/download/3PTB.pdb -o pdbs/3PTB.pdb
  ```
