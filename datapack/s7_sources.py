import json, glob, os
from common import *
L=['# SOURCES: NCASI Northern Forest manuscript data pack','',
   'Every output file in this folder, the script that wrote it (with md5) and its input files (with md5). Paths relative to the home folder on firebreather. '
   'Restricted inputs were read server side only; every output is a state-level aggregate. No plot coordinates, plot CNs or PLT_CN are written.','']
for f in sorted(glob.glob(f'{SRC}/*.json')):
    r=json.load(open(f)); out=f'{DATA}/{r["output"]}'
    L += [f'## {r["output"]}','',f'- output md5: {md5(out) if os.path.exists(out) else "MISSING"}',f'- script: {r["script"]} (md5 {r["script_md5"]})',f'- written: {r["written"]}']
    if r['notes']: L.append(f'- notes: {r["notes"]}')
    L.append('- inputs:')
    L += [f'  - {i["path"]} (md5 {i["md5"]})' for i in r['inputs']]
    L.append('')
L += ['## Helper scripts','',f'- {os.path.relpath(ROOT,HOME)}/tools/common.py (md5 {md5(ROOT+"/tools/common.py")})',
      f'- {os.path.relpath(ROOT,HOME)}/tools/run_queue.sh (md5 {md5(ROOT+"/tools/run_queue.sh")})','']
open(f'{DATA}/SOURCES.md','w').write('\n'.join(L))
status('step7 SOURCES.md written')
