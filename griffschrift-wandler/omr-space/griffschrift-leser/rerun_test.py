import sys,json;sys.path.insert(0,'/tmp/claude-0/gsr');from fuchsread import check
out={}
for k in json.load(open('/tmp/claude-0/gsr/midi_test.json')):
    r=check(f'/tmp/claude-0/fuchs/fuchs/{k}/');out[k]=r;print(k,r['n'],r['paired'],r['exact'],flush=True)
json.dump(out,open('/tmp/claude-0/gsr/midicheck_f.json','w'),default=lambda o:o.item() if hasattr(o,'item') else str(o))
