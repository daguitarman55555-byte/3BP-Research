# Fifth-stage evidence: the hidden finite cover for signed-area histories

Read `REPORT.md` first (it opens with the DECISIVE RESULT and the referee audit).

## Layout
* `code/` – all source code (Python 3.11; numpy, scipy, numba, python-flint, mpmath; msolve 0.6.5 optional)
  * `jets.py`, `jets_massvar.py` – exact area-jet polynomials on the lifted variety V (FLINT)
  * `numeval.py`, `track.py`, `track_robust.py`, `track_mass.py` – compiled evaluators / homotopy trackers
  * `monodromy.py`, `orbit.py`, `orbit_extend.py`, `orbit_extend_r.py`, `orbit_grow.py` – fibre construction by monodromy
  * `membership.py` – random-point completeness estimator
  * `cert_arb.py`, `run_arb.py`, `run_arb_new.py`, `analyze_cert.py` – rigorous ball-arithmetic Krawczyk certificates
  * `certify.py`, `run_certify.py` – earlier double-precision majorant test (superseded by arb, kept for record)
  * `exact_rank.py`, `check_engine.py`, `audit_equilateral.py`, `audit_37.py` – Phase-1 audits
  * `symmetry_tests.py`, `invariants.py`, `regimes.py`, `mass_transport*.py`, `move_fibre.py` – Phases 5,6,9,11
  * `make_msolve.py` – exact F_p Groebner attempt input generator
* `certificates/` – exact certificates and JSON summaries
* `logs/` – raw run logs
* `data_*` – jet polynomials (pickled exact rationals), fibre point sets (npz)

## Reproduce (approximate wall times on 4 cores)
```
cd code
python3 jets.py 1 2 3 8                      # seconds
python3 orbit.py 123 7 6 1.0                 # ~45 min, false plateau 19893
python3 orbit_extend.py ../data_orbit_123_s7.npz 101 4 1.0 ../data_orbit_123_ext1.npz
python3 orbit_extend.py ../data_orbit_123_ext1.npz 202 4 1.5 ../data_orbit_123_ext2.npz
python3 orbit_extend_r.py ../data_orbit_123_ext2.npz 303 4 1.0 ../data_orbit_123_ext3.npz
python3 membership.py ../data_orbit_123_ext3.npz 2000 2 ../certificates/membership_ext3_s2.json
python3 orbit_grow.py ../data_orbit_123_ext3.npz ../certificates/membership_ext3_s2.json,../certificates/membership_ext2_s1.json 404 6 ../data_orbit_123_ext4.npz
python3 run_arb.py ../data_orbit_123_s7.npz ../certificates/arb_fibre_123.pkl          # ~25 min
python3 run_arb_new.py ../data_orbit_123_ext4.npz ../certificates/arb_fibre_123.pkl ../certificates/arb_fibre_123_final.pkl
python3 analyze_cert.py ../certificates/arb_fibre_123_final.pkl ../certificates/fibre_123_final_summary.json
sh run_batch.sh                              # membership, regimes, audits, mass transports
```
Monodromy runs are randomized; seeds are fixed in the commands, but thread scheduling can change
which (rare) paths fail, so intermediate counts may differ slightly. The certificate is
re-checkable independently: `certificates/arb_fibre_123_final.pkl` contains every centre, and
`cert_arb.krawczyk` re-proves each box from the exact polynomials.
