import hashlib, json, os, datetime
HOME=os.path.expanduser('~')
ROOT=f'{HOME}/perseus_me/tracks/ncasi_ms'
DATA=f'{ROOT}/data'
SRC=f'{DATA}/_sources'
STATES=['ME','NH','VT','MI','WI','MN']
os.makedirs(SRC,exist_ok=True)
def md5(p):
    h=hashlib.md5()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()
def record(output, inputs, script, notes=''):
    rec={'output':output,'script':os.path.relpath(script,HOME) if script else '',
         'script_md5':md5(script) if script else '',
         'inputs':[{'path':os.path.relpath(p,HOME),'md5':md5(p)} for p in inputs],
         'notes':notes,'written':datetime.datetime.now().isoformat(timespec='seconds')}
    with open(f'{SRC}/{output}.json','w') as f: json.dump(rec,f,indent=1)
def status(msg):
    with open(f'{ROOT}/STATUS.md','a') as f:
        f.write(f'- {datetime.datetime.now().isoformat(timespec="minutes")} {msg}\n')
