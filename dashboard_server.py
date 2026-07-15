#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dashboard_server.py
===================
Serves the MOSSHOOS live dashboard from the server, with a real API that
RE-RUNS the models on demand and lets you download every input and output.

Standard library only -- no pip install required.

Usage
-----
    /usr/bin/python3 dashboard_server.py --port 8050
    then open  http://localhost:8050  (Cursor/VS Code forwards the port
    automatically; otherwise use an ssh tunnel:
      ssh -L 8050:localhost:8050 steeven@srv-ortad )

Routes
------
  GET  /                       dashboard
  GET  /api/results            current results JSON
  GET  /api/sweep              storage-sensitivity curve JSON
  POST /api/run                re-run the pipeline with supplied parameters
  GET  /api/download/<name>    inflows_omvs.csv | results.json | results.csv |
                               sweep.json | plan4res_storage.csv |
                               plan4res_bellman.csv | plan4res_inflow.csv
"""
import argparse
import json
import os
import subprocess
import sys
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse, parse_qs

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, 'results', 'comparison_results.json')
SWEEP = os.path.join(HERE, 'results', 'storage_sweep.json')
HTML = os.path.join(HERE, 'dashboard_app.html')
PY = os.environ.get('PYTHON', '/usr/bin/python3')
P4R_DEFAULT = os.path.expanduser('~/plan4res_work/data/toyDataset')

# run state shared with the browser
STATE = {'running': False, 'log': [], 'started': None, 'finished': None,
         'rc': None, 'params': None}
LOCK = threading.Lock()


def _p4r_file(*parts):
    return os.path.join(P4R_DEFAULT, *parts)


DOWNLOADS = {
    'inflows_omvs.csv':      (os.path.join(HERE, 'data', 'OMVS_inflows.csv'), 'text/csv'),
    'results.json':          (RESULTS, 'application/json'),
    'results.csv':           (os.path.join(HERE, 'results', 'comparison_results.csv'), 'text/csv'),
    'sweep.json':            (SWEEP, 'application/json'),
    'plan4res_storage.csv':  (_p4r_file('csv_simul', 'SS_SeasonalStorage.csv'), 'text/csv'),
    'plan4res_bellman.csv':  (_p4r_file('csv_simul', 'bellmanvalues.csv'), 'text/csv'),
    'plan4res_zonevalues.csv': (_p4r_file('csv_simul', 'ZV_ZoneValues.csv'), 'text/csv'),
    'plan4res_inflow.csv':   (_p4r_file('TimeSeries',
                              'AggregatedTimeSerie_Hydro-Run of River_RoGonDor.csv'), 'text/csv'),
}


def run_pipeline(params):
    """Run compare_plan4res.py with the given params (blocking)."""
    with LOCK:
        STATE.update(running=True, log=[], started=time.time(),
                     finished=None, rc=None, params=params)

    cmd = [PY, os.path.join(HERE, 'compare_plan4res.py')]
    if params.get('source') == 'plan4res':
        cmd += ['--dataset', params.get('dataset') or P4R_DEFAULT,
                '--reservoir', params.get('reservoir', 'RoGonDor')]
    else:
        cmd += ['--inflow', os.path.join(HERE, 'data', 'OMVS_inflows.csv')]
    cmd += ['--n_iter', str(int(params.get('n_iter', 150))),
            '--n_mc', str(int(params.get('n_mc', 200))),
            '--n_bootstrap', str(int(params.get('n_bootstrap', 150))),
            '--n_clusters', str(int(params.get('n_clusters', 10)))]

    with LOCK:
        STATE['log'].append('$ ' + ' '.join(cmd))
    try:
        proc = subprocess.Popen(cmd, cwd=HERE, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True, bufsize=1)
        for line in proc.stdout:
            line = line.rstrip('\n')
            with LOCK:
                STATE['log'].append(line)
                if len(STATE['log']) > 800:
                    STATE['log'] = STATE['log'][-800:]
        proc.wait()
        rc = proc.returncode
    except Exception as e:
        with LOCK:
            STATE['log'].append(f'ERROR: {e}')
        rc = 1
    with LOCK:
        STATE.update(running=False, finished=time.time(), rc=rc)


class Handler(BaseHTTPRequestHandler):
    def log_message(self, fmt, *a):
        pass  # quiet

    def _send(self, code, ctype, body, extra=None):
        if isinstance(body, str):
            body = body.encode('utf-8')
        self.send_response(code)
        self.send_header('Content-Type', ctype)
        self.send_header('Content-Length', str(len(body)))
        self.send_header('Cache-Control', 'no-store')
        for k, v in (extra or {}).items():
            self.send_header(k, v)
        self.end_headers()
        self.wfile.write(body)

    def _json(self, obj, code=200):
        self._send(code, 'application/json', json.dumps(obj))

    def do_GET(self):
        u = urlparse(self.path)
        p = u.path

        if p in ('/', '/index.html'):
            if not os.path.exists(HTML):
                return self._send(500, 'text/plain',
                                  f'dashboard_app.html not found at {HTML}')
            with open(HTML, 'r') as f:
                return self._send(200, 'text/html; charset=utf-8', f.read())

        if p == '/api/results':
            if not os.path.exists(RESULTS):
                return self._json({'error': 'no results yet - run the pipeline'}, 404)
            with open(RESULTS) as f:
                return self._send(200, 'application/json', f.read())

        if p == '/api/sweep':
            if not os.path.exists(SWEEP):
                return self._json({'error': 'no sweep yet'}, 404)
            with open(SWEEP) as f:
                return self._send(200, 'application/json', f.read())

        if p == '/api/status':
            with LOCK:
                s = dict(running=STATE['running'], rc=STATE['rc'],
                         params=STATE['params'],
                         elapsed=(round((STATE['finished'] or time.time())
                                        - STATE['started'], 1)
                                  if STATE['started'] else None),
                         log=STATE['log'][-200:])
            return self._json(s)

        if p == '/api/manifest':
            out = []
            for name, (path, _ct) in DOWNLOADS.items():
                out.append(dict(name=name, exists=os.path.exists(path),
                                size=(os.path.getsize(path)
                                      if os.path.exists(path) else 0),
                                path=path))
            return self._json(out)

        if p.startswith('/api/download/'):
            name = p.split('/api/download/', 1)[1]
            if name not in DOWNLOADS:
                return self._json({'error': 'unknown file'}, 404)
            path, ctype = DOWNLOADS[name]
            if not os.path.exists(path):
                return self._json({'error': f'missing: {path}'}, 404)
            with open(path, 'rb') as f:
                data = f.read()
            return self._send(200, ctype, data,
                              {'Content-Disposition': f'attachment; filename="{name}"'})

        return self._json({'error': 'not found'}, 404)

    def do_POST(self):
        u = urlparse(self.path)
        if u.path != '/api/run':
            return self._json({'error': 'not found'}, 404)
        with LOCK:
            if STATE['running']:
                return self._json({'error': 'a run is already in progress'}, 409)
        n = int(self.headers.get('Content-Length', 0))
        try:
            params = json.loads(self.rfile.read(n) or b'{}')
        except Exception:
            params = {}
        threading.Thread(target=run_pipeline, args=(params,), daemon=True).start()
        return self._json({'started': True, 'params': params})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--port', type=int, default=8050)
    ap.add_argument('--host', default='127.0.0.1')
    args = ap.parse_args()
    srv = ThreadingHTTPServer((args.host, args.port), Handler)
    print(f'MOSSHOOS dashboard serving on http://{args.host}:{args.port}')
    print(f'  results : {RESULTS}')
    print(f'  html    : {HTML}')
    print(f'  python  : {PY}')
    print('Ctrl+C to stop.')
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        print('\nstopped.')


if __name__ == '__main__':
    main()
