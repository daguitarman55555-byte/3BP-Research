while ps -eo args | grep -q "^sh run_batch.sh"; do sleep 20; done
python3 regimes_b.py ../data_orbit_123_ext4.npz ../certificates/regimes_123_b.json > ../logs/regimes_b.out 2>&1
python3 audit_37b.py > ../logs/audit_37b.out 2>&1
