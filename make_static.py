#!/usr/bin/env python3
"""
make_static.py
==============
Bake the live dashboard into a STANDALONE dashboard.html: all results embedded,
downloads become inline data: URIs, no server, no Python. Send it to anyone.

    python3 make_static.py                  # -> dashboard.html
    python3 make_static.py -o report.html
"""
import argparse, base64, json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))

FILES = {
    'inflows_omvs.csv':        ('data/OMVS_inflows.csv', 'text/csv'),
    'results.json':            ('results/comparison_results.json', 'application/json'),
    'results.csv':             ('results/comparison_results.csv', 'text/csv'),
    'sweep.json':              ('results/storage_sweep.json', 'application/json'),
    'plan4res_storage.csv':    (os.path.expanduser(
        '~/plan4res_work/data/toyDataset/csv_simul/SS_SeasonalStorage.csv'), 'text/csv'),
    'plan4res_bellman.csv':    (os.path.expanduser(
        '~/plan4res_work/data/toyDataset/csv_simul/bellmanvalues.csv'), 'text/csv'),
    'plan4res_zonevalues.csv': (os.path.expanduser(
        '~/plan4res_work/data/toyDataset/csv_simul/ZV_ZoneValues.csv'), 'text/csv'),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-o', '--output', default=os.path.join(HERE, 'dashboard.html'))
    ap.add_argument('--app', default=os.path.join(HERE, 'dashboard_app.html'))
    args = ap.parse_args()

    if not os.path.exists(args.app):
        sys.exit(f'not found: {args.app}')
    html = open(args.app).read()

    rp = os.path.join(HERE, 'results', 'comparison_results.json')
    if not os.path.exists(rp):
        sys.exit('no results yet - run the pipeline first')
    R = json.load(open(rp))
    sp = os.path.join(HERE, 'results', 'storage_sweep.json')
    SW = json.load(open(sp)) if os.path.exists(sp) else None

    # embed every available file as a data: URI
    man, embedded = [], {}
    for name, (path, ctype) in FILES.items():
        if not os.path.isabs(path):
            path = os.path.join(HERE, path)
        if os.path.exists(path):
            with open(path, 'rb') as f:
                raw = f.read()
            embedded[name] = (f'data:{ctype};base64,'
                              + base64.b64encode(raw).decode())
            man.append(dict(name=name, exists=True, size=len(raw)))
        else:
            man.append(dict(name=name, exists=False, size=0))

    # replace the fetch-based loader with baked-in data
    stub = f"""async function loadAll(){{
  R = {json.dumps(R)};
  SW = {json.dumps(SW)};
  render(); if(SW) renderSweep(); renderDownloads({json.dumps(man)});
}}
const EMBEDDED = {json.dumps(embedded)};
"""
    html2, n = re.subn(r'async function loadAll\(\)\{.*?\n\}\n', stub, html,
                       count=1, flags=re.S)
    if not n:
        sys.exit('could not patch loadAll() - dashboard_app.html changed?')
    html = html2

    # downloads point at the embedded data instead of the API
    html = html.replace('href="/api/download/${f.name}"',
                        'href="${EMBEDDED[f.name]||\'#\'}" download="${f.name}"')

    # the Run tab needs a server; say so plainly rather than showing dead controls
    html = html.replace('<div id="runAlert"></div>',
        '<div class="alert aw"><b>Static export.</b> This file has the results baked in and '
        'works offline with no server. To change parameters and re-run the models, run '
        '<code>python3 dashboard_server.py --port 8050</code> on a machine with the '
        'repository.</div><div id="runAlert"></div>', 1)
    html = html.replace("$('#btn-run').onclick=async()=>{",
                        "$('#btn-run').disabled=true;\n$('#btn-refresh').style.display='none';\n"
                        "$('#btn-run').onclick=async()=>{return;\n", 1)

    with open(args.output, 'w') as f:
        f.write(html)
    kb = os.path.getsize(args.output) / 1024
    print(f'Static dashboard: {args.output}  ({kb:.0f} KB)')
    print(f'  results embedded : {R["source"]}')
    print(f'  files embedded   : {sum(1 for m in man if m["exists"])}/{len(man)}')
    print('  opens in any browser, offline, no server needed')


if __name__ == '__main__':
    main()
