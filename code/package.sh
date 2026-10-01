cd /home/user/3BP-Research
python3 - <<'PY'
import os,hashlib
rows=[]
for root,_,fs in os.walk('.'):
    if '.git' in root or 'evidence_zip' in root: continue
    for f in sorted(fs):
        p=os.path.join(root,f)[2:]
        if p in ('SHA256SUMS.txt','MANIFEST.md') or p.endswith('.zip') or '__pycache__' in p: continue
        rows.append((p,os.path.getsize(p),hashlib.sha256(open(p,'rb').read()).hexdigest()))
rows.sort()
open('SHA256SUMS.txt','w').write(''.join(f'{h}  {p}\n' for p,s,h in rows))
desc={'code/':'source code','certificates/':'exact certificates / JSON summaries','logs/':'raw logs','data_':'data (exact jet polynomials, fibre point sets)'}
with open('MANIFEST.md','w') as m:
    m.write('# MANIFEST\n\n| file | bytes | role |\n|---|---|---|\n')
    for p,s,h in rows:
        role=next((v for k,v in desc.items() if p.startswith(k)),'top-level document')
        m.write(f'| `{p}` | {s} | {role} |\n')
PY
rm -f Fifth_Stage_Area_History_Evidence.zip
zip -qr Fifth_Stage_Area_History_Evidence.zip . -x '.git/*' '*__pycache__*' '*.zip'
ls -la Fifth_Stage_Area_History_Evidence.zip
