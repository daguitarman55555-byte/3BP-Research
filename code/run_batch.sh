cd /home/user/3BP-Research/code
python3 regimes.py ../data_orbit_123_ext4.npz ../certificates/regimes_123.json > ../logs/regimes.out 2>&1
sed -i 's#data_orbit_123_s7.npz#data_orbit_123_ext4.npz#' audit_37.py
python3 audit_37.py > ../logs/audit_37.out 2>&1
python3 mass_transport.py ../data_orbit_123_ext4.npz > ../logs/mass_transport.out 2>&1
python3 invariants.py ../data_orbit_123_ext4.npz ../certificates/invariant_degrees_final.json > ../logs/invariants.out 2>&1
echo ALLDONE >> ../logs/regimes.out
python3 mass_transport_122.py ../data_orbit_123_ext4.npz > ../logs/mass_122.out 2>&1
