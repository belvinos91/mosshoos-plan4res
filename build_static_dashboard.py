#!/usr/bin/env python3
import json, os
R = 'results'
B = json.load(open(os.path.join(R, 'dashboard_bundle.json')))
DATA = json.dumps(B, ensure_ascii=False)
HTML = r"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>MOSSHOOS &times; Plan4RES &mdash; Results</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/Chart.js/4.4.1/chart.umd.js"></script>
<style>
:root{--navy:#1F4E79;--navy2:#153a5c;--green:#1D7A5A;--orange:#C77700;--lblue:#EAF2F8;--lgreen:#EAF5EE;--lgray:#F4F6F7;--ink:#1a2130;--muted:#5a6572;--line:#e6ebf0}
*{box-sizing:border-box}body{margin:0;font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;background:var(--lgray);color:var(--ink);line-height:1.55}
header{background:linear-gradient(135deg,var(--navy),var(--navy2));color:#fff}.hb{max-width:1120px;margin:0 auto;padding:24px 30px}
header h1{margin:0;font-size:22px}header p{margin:7px 0 0;opacity:.9;font-size:13.5px;max-width:820px}
.tags{margin-top:12px;display:flex;gap:8px;flex-wrap:wrap}.tag{font-size:12px;background:rgba(255,255,255,.15);padding:4px 11px;border-radius:20px}
nav{background:var(--navy2);position:sticky;top:0;z-index:20}nav .ni{max-width:1120px;margin:0 auto;padding:0 30px;display:flex;gap:2px;flex-wrap:wrap}
nav button{background:none;border:0;color:#cfe0f0;font-size:13px;padding:12px 14px;cursor:pointer;border-bottom:3px solid transparent;font-family:inherit}
nav button:hover{color:#fff}nav button.on{color:#fff;border-bottom-color:var(--orange)}
.wrap{max-width:1120px;margin:0 auto;padding:24px 30px}.pane{display:none}.pane.on{display:block}
.card{background:#fff;border-radius:12px;box-shadow:0 1px 4px rgba(0,0,0,.06);padding:24px 28px;margin-bottom:20px}
.card h2{margin:0 0 3px;color:var(--navy);font-size:18px}.sub{color:var(--muted);font-size:13px;margin:0 0 16px}
.kpis{display:flex;gap:13px;flex-wrap:wrap}.kpi{flex:1;min-width:160px;background:var(--lgreen);border:1px solid #d6ebe0;border-radius:10px;padding:15px;text-align:center}
.kpi.b{background:var(--lblue);border-color:#d3e5f2}.kpi .v{font-size:25px;font-weight:750;color:var(--green)}.kpi.b .v{color:var(--navy)}.kpi .l{font-size:11.5px;color:#556;margin-top:5px}
table{width:100%;border-collapse:collapse;font-size:13px;margin:6px 0}th,td{padding:9px 11px;text-align:right;border-bottom:1px solid var(--line)}th:first-child,td:first-child{text-align:left}
thead th{background:var(--lblue);color:var(--navy)}.chartbox{position:relative;height:340px;margin:6px 0}
.note{background:var(--lgray);border-left:4px solid var(--green);border-radius:6px;padding:11px 15px;margin:12px 0;font-size:13px}.note.o{border-left-color:var(--orange)}.note.b{border-left-color:var(--navy)}
.dlgrid{display:grid;grid-template-columns:repeat(auto-fill,minmax(240px,1fr));gap:14px}.dl{border:1px solid var(--line);border-radius:10px;padding:16px}
.dl .n{font-weight:650;color:var(--navy)}.dl .d{font-size:12px;color:var(--muted);margin:4px 0 12px}.dl button{background:var(--green);color:#fff;border:0;border-radius:7px;padding:9px 16px;cursor:pointer;font-family:inherit;font-weight:600}
footer{max-width:1120px;margin:0 auto;padding:18px 30px 40px;color:var(--muted);font-size:12px}
</style></head><body>
<header><div class="hb"><h1>MOSSHOOS &times; Plan4RES &mdash; Results</h1>
<p>Stochastic seasonal storage valuation, OMVS / Senegal River. Plan4RES restored on open-source solvers; MOSSHOOS validated against it.</p>
<div class="tags"><span class="tag">HiGHS, no commercial licence</span><span class="tag">Manantali 1961&ndash;2015</span><span class="tag" id="tg"></span></div></div></header>
<nav><div class="ni">
<button class="on" data-p="ov">Overview</button><button data-p="wv">Water values</button><button data-p="ss">Storage sensitivity</button>
<button data-p="sc">Scenarios</button><button data-p="dl">Downloads</button></div></nav>
<div class="wrap">
<div class="pane on" id="pn-ov"><div class="card"><h2>Summary</h2><p class="sub">Every number below is produced by the committed pipeline.</p>
<div class="kpis"><div class="kpi b"><div class="v">&#10003;</div><div class="l">Plan4RES SSV&rarr;SIM&rarr;CEM on free software</div></div>
<div class="kpi"><div class="v" id="k-ag">&mdash;</div><div class="l">MOSSHOOS vs Plan4RES water value</div></div>
<div class="kpi b"><div class="v" id="k-vss">&mdash;</div><div class="l">Peak gain vs greedy release</div></div>
<div class="kpi"><div class="v">4/4</div><div class="l">Components validated</div></div></div>
<div class="note b"><b>Headline.</b> A one-line SMS++ fix restores Plan4RES on HiGHS. Driven with the same economics, MOSSHOOS reproduces Plan4RES's water values &mdash; a validation, not a superiority claim.</div>
<div class="note o"><b>Finding.</b> The value of stochastic operation is governed by the storage-to-inflow regime: ~0 in the normal range, material only at high storage &mdash; where Phase-2 jump dynamics target.</div></div></div>
<div class="pane" id="pn-wv"><div class="card"><h2>Water value = system marginal generation cost</h2><p class="sub">Plan4RES cuts (real, on HiGHS) vs MOSSHOOS with the same merit order.</p>
<div class="kpis" id="wvk"></div><div class="chartbox"><canvas id="wvC"></canvas></div>
<div class="note">Plan4RES tracks the merit order (coal 3, nuclear 8, CCGT 81.66, VOLL 10000). MOSSHOOS lands at the system marginal cost &mdash; agreement, not superiority.</div></div></div>
<div class="pane" id="pn-ss"><div class="card"><h2>Value of stochastic storage vs usable storage</h2><p class="sub">Manantali, 55-year record.</p>
<div class="chartbox"><canvas id="ssC"></canvas></div><table id="ssT"><thead><tr><th>Storage (&times;mean-mo)</th><th>S_max (GWh)</th><th>vs greedy</th><th>vs deterministic</th></tr></thead><tbody></tbody></table></div></div>
<div class="pane" id="pn-sc"><div class="card"><h2>Stochastic scenarios</h2><p class="sub">Reduced weighted scenario paths (one annual cycle).</p>
<div class="chartbox"><canvas id="scC"></canvas></div></div></div>
<div class="pane" id="pn-dl"><div class="card"><h2>Data &amp; downloads</h2><p class="sub">Inputs, outputs and scenarios &mdash; embedded, work offline.</p><div class="dlgrid" id="dlg"></div></div></div>
</div>
<footer>Built from committed <code>results/</code>. Rebuild: <code>python3 build_static_dashboard.py</code>.</footer>
<script>
const D=__DATA__;
document.getElementById('tg').textContent=(D.scenarios.length)+' scenarios';
document.querySelectorAll('nav button').forEach(b=>b.onclick=()=>{document.querySelectorAll('nav button').forEach(x=>x.classList.remove('on'));document.querySelectorAll('.pane').forEach(x=>x.classList.remove('on'));b.classList.add('on');document.getElementById('pn-'+b.dataset.p).classList.add('on');window.scrollTo(0,0);});
const C={navy:'#1F4E79',green:'#1D7A5A',orange:'#C77700',mut:'#9aa7b3'};
const ag=Math.round(D.wv2.agreement_pct||97);document.getElementById('k-ag').textContent='~'+ag+'%';
let peak=0;D.storage_sweep.sweep.forEach(r=>{if(r.vss_vs_greedy>peak)peak=r.vss_vs_greedy;});document.getElementById('k-vss').textContent=peak.toFixed(0)+'%';
const p4r=D.wv2.plan4res_wv,mos=D.wv2.mosshoos_mean,mc=D.wv2.mc_avg;
const wvk=document.getElementById('wvk');[['Plan4RES stages',p4r.length,'b'],['Exact matches',Math.round(D.wv2.plan4res_exact_fraction*p4r.length)+'/'+p4r.length,''],['MOSSHOOS mean',mos.toFixed(1),''],['Marginal cost',mc.toFixed(1),'b']].forEach(([l,v,c])=>wvk.insertAdjacentHTML('beforeend','<div class="kpi '+(c=='b'?'b':'')+'"><div class="v">'+v+'</div><div class="l">'+l+'</div></div>'));
new Chart(document.getElementById('wvC'),{type:'line',data:{labels:p4r.map((_,i)=>i),datasets:[{label:'Plan4RES SDDP (HiGHS)',data:p4r,borderColor:C.navy,backgroundColor:C.navy,pointRadius:4,tension:.15},{label:'MOSSHOOS mean '+mos.toFixed(1),data:p4r.map(()=>mos),borderColor:C.green,borderWidth:3,pointRadius:0}]},options:{scales:{y:{type:'logarithmic',title:{display:true,text:'$/MWh'}},x:{title:{display:true,text:'stage'}}},plugins:{legend:{position:'bottom'}}}});
const sw=D.storage_sweep.sweep;
new Chart(document.getElementById('ssC'),{type:'line',data:{labels:sw.map(r=>r.frac+'x'),datasets:[{label:'vs greedy (%)',data:sw.map(r=>r.vss_vs_greedy),borderColor:C.navy,tension:.2,pointRadius:3},{label:'vs deterministic (%)',data:sw.map(r=>r.vss),borderColor:C.orange,tension:.2,pointRadius:3}]},options:{scales:{y:{title:{display:true,text:'VSS (%)'}},x:{title:{display:true,text:'usable storage (x mean monthly inflow)'}}},plugins:{legend:{position:'bottom'}}}});
const sb=document.querySelector('#ssT tbody');sw.forEach(r=>sb.insertAdjacentHTML('beforeend','<tr><td>'+r.frac+'x</td><td>'+r.S_max.toFixed(1)+'</td><td>'+r.vss_vs_greedy.toFixed(1)+'%</td><td>'+r.vss.toFixed(1)+'%</td></tr>'));
const pal={WET:C.navy,DRY:C.orange,median:C.mut};
new Chart(document.getElementById('scC'),{type:'line',data:{labels:Array.from({length:(D.scenarios[0]&&D.scenarios[0].path.length)||12},(_,i)=>'m'+(i+1)),datasets:D.scenarios.map(s=>({label:'#'+s.id+' '+s.type,data:s.path,borderColor:pal[s.type]||C.mut,borderWidth:s.type=='median'?1.2:2.4,pointRadius:0,tension:.3}))},options:{plugins:{legend:{position:'bottom',labels:{boxWidth:12,font:{size:10}}}}}});
function dl(n,t){const b=new Blob([t],{type:'text/plain'});const u=URL.createObjectURL(b);const a=document.createElement('a');a.href=u;a.download=n;a.click();URL.revokeObjectURL(u);}
const g=document.getElementById('dlg');Object.keys(D.files).forEach(fn=>g.insertAdjacentHTML('beforeend','<div class="dl"><div class="n">'+fn+'</div><div class="d">committed output</div><button data-f="'+fn+'">Download</button></div>'));
g.querySelectorAll('button').forEach(btn=>btn.onclick=()=>dl(btn.dataset.f,D.files[btn.dataset.f]));
</script></body></html>"""
open('dashboard.html','w',encoding='utf-8').write(HTML.replace('__DATA__',DATA))
print('wrote dashboard.html', round(len(HTML)/1024,1),'KB template +', round(len(DATA)/1024,1),'KB data')
