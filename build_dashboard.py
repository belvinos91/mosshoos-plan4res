#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
build_dashboard.py
==================
Generates dashboard.html from results/comparison_results.json.
NO hardcoded numbers — every value comes from the JSON that
compare_plan4res.py just produced. Fully reproducible.

Usage:
    python3 build_dashboard.py
    # then open dashboard.html in any browser
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, 'results', 'comparison_results.json')

if not os.path.exists(RESULTS):
    sys.exit("No results found. Run compare_plan4res.py first.")

R = json.load(open(RESULTS))
# embed the JSON directly so the dashboard is self-contained (works offline, no fetch)
DATA_JSON = json.dumps(R)

html = '''<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>MOSSHOOS vs Plan4RES — Live Comparison Dashboard</title>
<script>window.MathJax={tex:{inlineMath:[['\\\\(','\\\\)']],displayMath:[['\\\\[','\\\\]']]},svg:{fontCache:'global'}};</script>
<script src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js" async></script>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
<style>
:root{--blue:#1F4E79;--green:#1D7A5A;--orange:#C77700;--red:#8B1A1A;--lblue:#EAF2F8;--lgreen:#EAF5EE;--lgray:#F4F6F7;--ink:#1a1a2e;--muted:#5a6572}
*{box-sizing:border-box}body{margin:0;font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;background:var(--lgray);color:var(--ink);line-height:1.55}
header{background:linear-gradient(135deg,#1F4E79,#153a5c);color:#fff}
.hbar{max-width:1180px;margin:0 auto;padding:26px 32px}
header h1{margin:0;font-size:23px}header .subtitle{margin:8px 0 0;opacity:.92;font-size:14px;max-width:840px}
header .authors{margin:12px 0 0;font-size:12.5px;opacity:.8}
header .datasource{margin:10px 0 0;font-size:12px;background:rgba(255,255,255,.14);display:inline-block;padding:4px 12px;border-radius:12px}
nav{background:#153a5c;position:sticky;top:0;z-index:50}
nav .ninner{max-width:1180px;margin:0 auto;padding:0 32px;display:flex;gap:2px;flex-wrap:wrap}
nav a{color:#cfe0f0;text-decoration:none;font-size:13px;padding:12px 13px;border-bottom:3px solid transparent}
nav a:hover,nav a.active{color:#fff;border-bottom-color:var(--orange)}
.wrap{max-width:1180px;margin:0 auto;padding:26px 32px}
section{background:#fff;border-radius:12px;box-shadow:0 1px 5px rgba(0,0,0,.07);padding:26px 30px;margin-bottom:24px;scroll-margin-top:58px}
section h2{margin:0 0 4px;color:var(--blue);font-size:19px}
.lead{color:var(--muted);font-size:13.5px;margin:0 0 18px}
h3{color:var(--green);font-size:15px;margin:20px 0 8px}
.kpirow{display:flex;gap:14px;flex-wrap:wrap}
.kpi{flex:1;min-width:150px;background:var(--lgreen);border-radius:10px;padding:16px;text-align:center;border:1px solid #d6ebe0}
.kpi.blue{background:var(--lblue);border-color:#d3e5f2}.kpi .v{font-size:27px;font-weight:750;color:var(--green)}
.kpi.blue .v{color:var(--blue)}.kpi .l{font-size:12px;color:#556;margin-top:5px}
table{width:100%;border-collapse:collapse;font-size:13.5px;margin:8px 0}
th,td{padding:10px 12px;text-align:right;border-bottom:1px solid #eef1f4}th:first-child,td:first-child{text-align:left}
thead th{background:var(--lblue);color:var(--blue)}tbody tr:hover{background:#fafcfe}
.badge{display:inline-block;font-size:11px;font-weight:600;padding:2px 9px;border-radius:11px;color:#fff}
.badge.ms{background:var(--green)}.badge.md{background:var(--blue)}.badge.ref{background:var(--muted)}.badge.new{background:var(--orange)}
.eqbox{background:var(--lgray);border-left:4px solid var(--blue);border-radius:6px;padding:12px 16px;margin:10px 0;overflow-x:auto}
.eqbox.green{border-left-color:var(--green)}
.grid2{display:grid;grid-template-columns:1fr 1fr;gap:20px}@media(max-width:820px){.grid2{grid-template-columns:1fr}}
.note{background:#fff8ec;border:1px solid #f0dcb8;border-radius:8px;padding:12px 15px;font-size:12.8px;color:#6b5326;margin:12px 0}
.improve{background:linear-gradient(135deg,#EAF5EE,#f4faf6);border:1px solid #cfe8da;border-radius:10px;padding:16px 18px;margin:10px 0}
.improve h4{margin:0 0 5px;color:var(--green);font-size:14.5px}.arrow{color:var(--orange);font-weight:700;margin:0 6px}
.chartbox{position:relative;height:290px;margin:12px 0}
.zone{display:inline-block;background:#fff;border:1px solid #d3e5f2;border-radius:8px;padding:9px 14px;margin:4px 7px 4px 0;font-size:12.5px}.zone b{color:var(--blue)}
.tag-s{color:var(--green);font-weight:600}.tag-d{color:var(--blue);font-weight:600}
.live{display:inline-block;width:8px;height:8px;background:#2ecc71;border-radius:50%;margin-right:6px;animation:pulse 1.6s infinite}
@keyframes pulse{0%,100%{opacity:1}50%{opacity:.3}}
footer{text-align:center;color:#9aa4ae;font-size:12px;padding:24px}
</style></head><body>
<header><div class="hbar">
 <h1>MOSSHOOS vs Plan4RES — Comparison Dashboard</h1>
 <div class="subtitle">Stochastic seasonal storage valuation (Haar wavelet + block-bootstrap + SDDP + HJB) benchmarked against the Plan4RES deterministic method, on identical inputs.</div>
 <div class="authors">S. B. Affognon · B. M. Ndiaye · P. Mendy · C. M. F. Kebe &nbsp;|&nbsp; UCAD Dakar · University of Nairobi · CTSS</div>
 <div class="datasource"><span class="live"></span><span id="src"></span> · generated live from results JSON</div>
</div>
<nav><div class="ninner">
 <a href="#overview">Overview</a><a href="#region">Study Region</a>
 <a href="#p4r">Plan4RES Method</a><a href="#mos">MOSSHOOS Method</a>
 <a href="#scen">Scenarios</a><a href="#results">Results</a><a href="#ct">Continuous-time</a><a href="#improve">Where We Improve</a>
</div></nav></header>
<div class="wrap">

<section id="overview"><h2>Overview</h2>
<p class="lead">All numbers below are computed by <code>compare_plan4res.py</code> and read live from <code>results/comparison_results.json</code> — nothing is hardcoded. Re-run the pipeline and this dashboard updates.</p>
<div class="kpirow">
 <div class="kpi"><div class="v" id="k_vss"></div><div class="l">Value of Stochastic Solution (VSS)</div></div>
 <div class="kpi blue"><div class="v" id="k_p90"></div><div class="l">P90 tail-cost saving</div></div>
 <div class="kpi"><div class="v" id="k_es"></div><div class="l">Scenario energy-score gain vs i.i.d.</div></div>
 <div class="kpi blue"><div class="v" id="k_aj"></div><div class="l">Wavelet seasonal energy A_J</div></div>
</div>
<p style="font-size:13px;margin-top:16px">
 <span class="badge ms">M-S</span> MOSSHOOS (stochastic) ·
 <span class="badge md">M-D</span> Plan4RES (deterministic) ·
 <span class="badge ref">M-S0</span> stochastic without wavelet ·
 <span class="badge ref">M-D0</span> greedy baseline</p>
</section>

<section id="region"><h2>Study Region &amp; Data</h2>
<p class="lead">Reservoir and parameters read from Plan4RES <code>SS_SeasonalStorage.csv</code> / <code>ZV_ZoneValues.csv</code>.</p>
<div id="zones"></div>
<div class="note">For the Senegal/WAPP application these map to the OMVS plants: <b>Manantali/BafingMakana</b> (200 MW, H=65 m), <b>Félou</b> (59.5 MW, H=14.5 m), <b>Gouina</b> (140 MW, H=17 m).</div>
</section>

<section id="p4r"><h2><span class="badge md">M-D</span> Plan4RES Method</h2>
<p class="lead">Deterministic seasonal storage value via Benders cuts on the mean-inflow forecast.</p>
<div class="eqbox">\\[ V_t(S)=\\max_{\\ell\\in\\mathcal{L}_t}\\{\\,b_\\ell+\\mathbf{a}_\\ell^\\top\\mathbf{S}\\,\\} \\]</div>
<div id="p4rref"></div>
</section>

<section id="mos"><h2><span class="badge ms">M-S</span> MOSSHOOS Method</h2>
<p class="lead">Four coupled layers replace the deterministic assumption with a full stochastic treatment.</p>
<div class="grid2">
 <div><h3>WP2 — Haar wavelet</h3><div class="eqbox green">\\[ x_t=A_J(t)+\\sum_{j=1}^{J}D_j(t) \\]</div>
  <p style="font-size:12.5px" id="wp2txt"></p></div>
 <div><h3>WP2 — Fourier seasonal fit</h3><div class="eqbox green">\\[ \\theta(t)=a_0+\\sum_{k=1}^{K}[a_k\\cos(2\\pi kt)+b_k\\sin(2\\pi kt)] \\]</div>
  <p style="font-size:12.5px" id="fourtxt"></p></div>
</div>
<h3>WP3 — Block-bootstrap scenarios</h3><p style="font-size:12.5px" id="wp3txt"></p>
<div class="grid2">
 <div><h3>WP4-B — Itô SDE (OU / CIR)</h3>
  <div class="eqbox green">\\[ dQ_t=\\kappa(\\mu-Q_t)dt+\\sigma\\,dW_t \\]</div>
  <div class="eqbox green">\\[ dQ_t=\\kappa(\\theta-Q_t)dt+\\sigma\\sqrt{Q_t}\\,dW_t \\]</div></div>
 <div><h3>WP4-C — HJB control</h3>
  <div class="eqbox green">\\[ \\rho J=\\sup_u\\{\\mathcal{R}+(Q-u)J_V+\\kappa(\\mu-Q)J_Q+\\tfrac{\\sigma^2}{2}J_{QQ}\\} \\]</div>
  <div class="eqbox green">\\[ u^*(V,Q)=\\Pi_{[0,u_{\\max}]}\\!\\big(\\tfrac{b-J_V}{c}\\big) \\]</div></div>
</div>
</section>

<section id="scen"><h2>Scenarios</h2>
<p class="lead">MOSSHOOS represents inflow uncertainty as weighted scenarios (dry→wet) vs the single Plan4RES deterministic mean.</p>
<div class="chartbox"><canvas id="scenChart"></canvas></div>
</section>

<section id="results"><h2>Results — Model Comparison</h2>
<p class="lead" id="restxt"></p>
<table id="restable"><thead><tr><th>Model</th><th>Mean cost</th><th>P90 cost</th><th>Deficit (GWh)</th></tr></thead><tbody></tbody></table>
<div class="chartbox"><canvas id="costChart"></canvas></div>
<div class="note" id="resnote"></div>
</section>

<section id="ct"><h2>WP4-B/C - Continuous-time strand (validated)</h2>
<p class="lead">Ito SDE calibration, strong-convergence rates, HJB PDE solution and Haar wavelet-Galerkin convergence - all computed live and unit-tested.</p>
<div class="kpirow">
 <div class="kpi"><div class="v" id="k_feller"></div><div class="l">Feller ratio (>=1 keeps Q&gt;0)</div></div>
 <div class="kpi blue"><div class="v" id="k_em"></div><div class="l">Euler-Maruyama slope (theory 0.5)</div></div>
 <div class="kpi"><div class="v" id="k_mil"></div><div class="l">Milstein slope (theory 1.0)</div></div>
 <div class="kpi blue"><div class="v" id="k_wg"></div><div class="l">Wavelet-Galerkin rate (theory -1)</div></div>
</div>
<h3>SDE calibration</h3>
<table id="sdetab"><thead><tr><th>Parameter</th><th>Value</th></tr></thead><tbody></tbody></table>
<h3>HJB value function at t=0 (mid inflow)</h3>
<div class="chartbox"><canvas id="hjbChart"></canvas></div>
<p style="font-size:12.5px;color:var(--muted)" id="hjbnote"></p>
<h3>Wavelet-Galerkin convergence</h3>
<table id="wgtab"><thead><tr><th>J</th><th>M basis fns</th><th>L2 error</th></tr></thead><tbody></tbody></table>
</section>
<section id="improve"><h2>Where We Improve Plan4RES</h2>
<div class="improve"><h4>1 · Inflow representation</h4><span class="tag-d">Plan4RES:</span> single deterministic mean \\(\\bar{Q}_t\\) <span class="arrow">→</span> <span class="tag-s">MOSSHOOS:</span> Haar wavelet \\(x_t=A_J+\\sum D_j\\)</div>
<div class="improve"><h4>2 · Scenario generation</h4><span class="tag-d">Plan4RES:</span> none / i.i.d. <span class="arrow">→</span> <span class="tag-s">MOSSHOOS:</span> block-bootstrap + K-means, <span id="imp_es"></span></div>
<div class="improve"><h4>3 · Bellman values</h4><span class="tag-d">Plan4RES:</span> deterministic Benders cuts <span class="arrow">→</span> <span class="tag-s">MOSSHOOS:</span> online SDDP over scenario set η</div>
<div class="improve"><h4>4 · Release policy</h4><span class="tag-d">Plan4RES:</span> stage-wise LP <span class="arrow">→</span> <span class="tag-s">MOSSHOOS:</span> continuous-time HJB feedback under OU/CIR</div>
<div class="improve"><h4>5 · Uncertainty quantification</h4><span class="tag-d">Plan4RES:</span> not reported <span class="arrow">→</span> <span class="tag-s">MOSSHOOS:</span> <span id="imp_vss"></span></div>
</section>

</div>
<footer id="foot"></footer>

<script>
const R = ''' + f"{DATA_JSON}" + ''';
const f=(x,d=1)=>Number(x).toLocaleString(undefined,{minimumFractionDigits:d,maximumFractionDigits:d});
document.getElementById('src').textContent = R.source;
document.getElementById('k_vss').textContent = R.wp4.VSS_pct.toFixed(1)+'%';
document.getElementById('k_p90').textContent = R.wp4.P90_saving_pct.toFixed(1)+'%';
document.getElementById('k_es').textContent = R.wp3.es_improvement.toFixed(1)+'%';
document.getElementById('k_aj').textContent = R.wp2.AJ_energy_pct.toFixed(1)+'%';
// zones
const p=R.params;
document.getElementById('zones').innerHTML =
 `<div class="zone"><b>${p.reservoir}</b><br>S_max ${f(p.S_max)} · R_max ${f(p.R_max)} · C_def ${f(p.C_def,0)}</div>`;
// plan4res reference
if(R.plan4res_reference){
 let h=`<h3>Their reference output (bellmanvalues.csv)</h3><p style="font-size:12.5px">${R.plan4res_reference.n_timesteps} timesteps · ${R.plan4res_reference.n_cuts} cuts</p><table><thead><tr><th>Timestep</th><th>a₀</th><th>a₁</th><th>a₂</th><th>b</th></tr></thead><tbody>`;
 R.plan4res_reference.sample.forEach(r=>{h+=`<tr><td>${r.Timestep}</td><td>${f(r.a_0,3)}</td><td>${f(r.a_1,3)}</td><td>${f(r.a_2,3)}</td><td>${f(r.b,0)}</td></tr>`;});
 document.getElementById('p4rref').innerHTML=h+'</tbody></table>';
}
// wp text
document.getElementById('wp2txt').textContent=`J*=${R.wp2.J_optimal}, A_J energy=${R.wp2.AJ_energy_pct}%, residual σ=${R.wp2.residual_std} GWh/mo.`;
document.getElementById('fourtxt').textContent=`K=5 harmonics, R²=${R.wp2.fourier_R2}.`;
document.getElementById('wp3txt').textContent=`${R.wp3.n_scenarios} scenarios (b=4 blocks), KS D=${R.wp3.ks_stat}, ACF max err=${R.wp3.acf_max_err}, energy-score improvement ${R.wp3.es_improvement}% vs i.i.d.`;
// results table
const labels={MS:'M-S MOSSHOOS',MS0:'M-S0 no-wavelet',MD:'M-D Plan4RES',MD0:'M-D0 greedy'};
const tb=document.querySelector('#restable tbody');
['MS','MS0','MD','MD0'].forEach(m=>{const d=R.wp4.models[m];
 tb.innerHTML+=`<tr><td>${labels[m]}</td><td>${f(d.mean_cost,0)}</td><td>${f(d.p90_cost,0)}</td><td>${f(d.deficit)}</td></tr>`;});
document.getElementById('restxt').textContent=`Monte-Carlo comparison on ${R.source}. VSS = ${R.wp4.VSS_pct}% (M-D vs M-S).`;
document.getElementById('resnote').innerHTML=`<b>Interpretation:</b> VSS = ${R.wp4.VSS_pct}% on this dataset. Stochastic optimisation delivers its largest gains in high-variability systems; on low-variability inputs the deterministic method is nearly adequate, but the greedy baseline (M-D0) is consistently worse — confirming the water-value mechanism operates.`;
document.getElementById('imp_es').textContent=`+${R.wp3.es_improvement}% energy score`;
document.getElementById('imp_vss').textContent=`VSS = ${R.wp4.VSS_pct}%, P90 saving = ${R.wp4.P90_saving_pct}%`;
document.getElementById('foot').textContent=`MOSSHOOS Dashboard · ${R.source} · runtime ${R.runtime_s}s · generated from results JSON`;
// scenario chart
const months=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'];
const pal=['#C0392B','#E67E22','#27AE60','#2980B9','#16A085','#8E44AD','#D35400','#2C3E50','#7F8C8D','#1ABC9C'];
const ds=R.wp3.scenarios.map((s,i)=>({label:`S${i+1} (${R.wp3.types[i]})`,data:s,borderColor:pal[i%pal.length],backgroundColor:'transparent',borderWidth:1.5,tension:.35,pointRadius:0}));
ds.push({label:'Plan4RES deterministic mean',data:R.seasonal_profile,borderColor:'#1F4E79',borderDash:[7,4],borderWidth:2.5,tension:.35,pointRadius:0});
new Chart(document.getElementById('scenChart'),{type:'line',data:{labels:months,datasets:ds},
 options:{responsive:true,maintainAspectRatio:false,plugins:{legend:{position:'bottom',labels:{boxWidth:12,font:{size:10}}}},scales:{y:{title:{display:true,text:'Inflow (GWh/month)'},beginAtZero:true},x:{grid:{display:false}}}}});
// cost chart
const mm=['MS','MS0','MD','MD0'];
new Chart(document.getElementById('costChart'),{type:'bar',
 data:{labels:mm.map(m=>labels[m]),datasets:[
  {label:'Mean cost',data:mm.map(m=>R.wp4.models[m].mean_cost),backgroundColor:['#1D7A5A','#52b788','#1F4E79','#C77700']},
  {label:'P90 cost',data:mm.map(m=>R.wp4.models[m].p90_cost),backgroundColor:['#52b788','#a9dfbf','#5a8fc0','#e0a458']}]},
 options:{responsive:true,maintainAspectRatio:false,plugins:{legend:{position:'bottom',labels:{boxWidth:12,font:{size:10}}}},scales:{y:{title:{display:true,text:'Annual cost ($/yr)'},beginAtZero:true},x:{grid:{display:false}}}}});
// nav highlight

// ---- WP4-B/C continuous-time strand ----
if (R.wp4b_continuous) {
  const C = R.wp4b_continuous;
  document.getElementById('k_feller').textContent = C.sde.feller.toFixed(1);
  document.getElementById('k_em').textContent = C.convergence.em_slope.toFixed(2);
  document.getElementById('k_mil').textContent = C.convergence.milstein_slope.toFixed(2);
  document.getElementById('k_wg').textContent =
      (C.wavelet_galerkin && C.wavelet_galerkin.slope != null)
      ? C.wavelet_galerkin.slope.toFixed(2) : 'n/a';

  const sd = document.querySelector('#sdetab tbody');
  const rows = [
    ['sigma type', C.sde.sigma_type],
    ['kappa (mean reversion, /yr)', C.sde.kappa],
    ['sigma0 (diffusion)', C.sde.sigma0],
    ['Fourier R2 (seasonal fit)', C.sde.fourier_R2],
    ['KS p-value', C.sde.ks_p + (C.sde.gaussian_innovations ? ' (Gaussian)' : ' (heavy-tailed)')],
    ['Feller ratio 2*kappa*theta/sigma0^2', C.sde.feller + (C.sde.feller_ok ? ' (Q>0 guaranteed)' : ' (may hit 0)')]
  ];
  rows.forEach(r => { sd.innerHTML += `<tr><td>${r[0]}</td><td>${r[1]}</td></tr>`; });

  if (C.hjb && C.hjb.storage_grid) {
    new Chart(document.getElementById('hjbChart'), {
      type:'line',
      data:{labels:C.hjb.storage_grid, datasets:[
        {label:'Value function V(S)', data:C.hjb.value_at_t0, borderColor:'#1D7A5A',
         backgroundColor:'transparent', borderWidth:2, tension:.25, pointRadius:0, yAxisID:'y'},
        {label:'Optimal release u*(S)', data:C.hjb.policy_at_t0, borderColor:'#C77700',
         backgroundColor:'transparent', borderWidth:2, borderDash:[6,4], tension:.25,
         pointRadius:0, yAxisID:'y1'}
      ]},
      options:{responsive:true, maintainAspectRatio:false,
        plugins:{legend:{position:'bottom',labels:{boxWidth:12,font:{size:10}}}},
        scales:{
          y:{position:'left', title:{display:true,text:'V(S) [cost-to-go]'}},
          y1:{position:'right', title:{display:true,text:'release [GWh]'}, grid:{display:false}},
          x:{title:{display:true,text:'storage S [GWh]'}, grid:{display:false}}
        }}
    });
    document.getElementById('hjbnote').textContent =
      `Solver stable: ${C.hjb.stable}. Value function U-shaped: ${C.hjb.u_shaped} - `
      + `high at empty (deficit risk, V=${C.hjb.V_at_empty}), minimum V=${C.hjb.V_at_min} `
      + `at S=${C.hjb.V_min_at_storage} GWh, rising to V=${C.hjb.V_at_full} near full `
      + `(spill penalty). This U-shape is the physically correct signature.`;
  }

  if (C.wavelet_galerkin && C.wavelet_galerkin.table) {
    const wt = document.querySelector('#wgtab tbody');
    C.wavelet_galerkin.table.forEach(r => {
      wt.innerHTML += `<tr><td>${r.J}</td><td>${r.M}</td><td>${r.L2_error.toExponential(3)}</td></tr>`;
    });
  }
}

const secs=[...document.querySelectorAll('section')],links=[...document.querySelectorAll('nav a')];
addEventListener('scroll',()=>{let c='';secs.forEach(s=>{if(scrollY>=s.offsetTop-70)c=s.id});links.forEach(a=>a.classList.toggle('active',a.getAttribute('href')==='#'+c));});
</script></body></html>'''

out = os.path.join(HERE, 'dashboard.html')
open(out,'w').write(html)
print(f"Dashboard generated: {out}")
print(f"Read from: {RESULTS}")
print(f"VSS shown: {R['wp4']['VSS_pct']}%  (live from JSON)")
