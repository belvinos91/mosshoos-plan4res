#!/usr/bin/env python3
"""Regenerate the report figures from committed results (server-side)."""
import json, os
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
NAVY='#1F4E79'; GREEN='#1D7A5A'; ORANGE='#C77700'; MUT='#9aa7b3'
OUT=os.path.join(os.path.dirname(os.path.abspath(__file__)),'figures')
os.makedirs(OUT,exist_ok=True)
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Fig: water-value comparison (from committed JSON)
d=json.load(open(os.path.join(ROOT,'results','water_value_comparison_v2.json')))
p4r=np.array(d['plan4res_wv']); mos=d['mosshoos_mean']
fig,ax=plt.subplots(figsize=(8.2,4.6))
for v,l in [(10000,'VOLL 10000'),(113.4,'OCGT 113.4'),(81.665,'CCGT 81.66'),(8,'Nuclear 8.00'),(3,'Coal 3.00')]:
    ax.axhline(v,ls='--',lw=.7,color='0.7'); ax.text(len(p4r)-.3,v*1.1,l,fontsize=7.5,color='0.5',ha='right')
ax.plot(range(len(p4r)),p4r,'o-',color=NAVY,lw=1.6,ms=6,label='Plan4RES SDDP (real cuts, HiGHS)')
ax.axhline(mos,color=GREEN,lw=2.4,label=f'MOSSHOOS (real economics), mean {mos:.1f} $/MWh')
ax.set_yscale('log'); ax.set_xlabel('stage'); ax.set_ylabel('marginal water value [$/MWh]')
ax.set_title('Water value = system marginal generation cost (RoGonDor)')
ax.legend(loc='center left',fontsize=8); plt.tight_layout()
plt.savefig(os.path.join(OUT,'water_value_comparison_v2.png'),dpi=150); plt.close()

# Fig: storage sweep (from committed JSON)
sw=json.load(open(os.path.join(ROOT,'storage_sweep_results.json')))['sweep']
fr=[r['frac'] for r in sw]; vg=[r['vss_vs_greedy'] for r in sw]; vd=[r['vss'] for r in sw]
fig,ax=plt.subplots(figsize=(7.6,4.2))
ax.plot(fr,vg,'o-',color=NAVY,lw=2,label='vs unmanaged (greedy) release')
ax.plot(fr,vd,'s-',color=ORANGE,lw=2,label='vs deterministic-equivalent policy')
ax.axvspan(min(fr),4.2,color=MUT,alpha=.10)
ax.set_xlabel('usable storage (x mean monthly inflow)'); ax.set_ylabel('value of policy (% cost saving)')
ax.set_title('Value of seasonal storage valuation vs storage size - Manantali 1961-2015')
ax.legend(); ax.grid(alpha=.25); plt.tight_layout()
plt.savefig(os.path.join(OUT,'fig_storage_sweep.png'),dpi=150); plt.close()

# Fig: wavelet ablation (from committed JSON)
w=json.load(open(os.path.join(ROOT,'results','wavelet_ablation.json')))
means=[w['mean_wavelet'],w['mean_seasonal_only'],w['mean_naive']]
labels=['Wavelet\n(A_J+block detail)','Seasonal-only\n(climatology)','Naive\n(i.i.d. months)']
fig,ax=plt.subplots(figsize=(7.0,4.0))
bars=ax.bar(labels,means,color=[GREEN,NAVY,MUT],width=.58)
for b,m in zip(bars,means): ax.text(b.get_x()+b.get_width()/2,m+250,f'{m:,.0f}',ha='center',fontsize=9)
lo,hi=w['ci']; ax.set_ylabel('mean out-of-sample cost ($/annual cycle)')
ax.set_title(f"Scenario ablation, 6x storage, {w['n_eval']} real years - naive-wavelet "
             f"{w['diff_naive_minus_wavelet']:+,.0f} [CI {lo:+,.0f},{hi:+,.0f}]")
ax.set_ylim(min(means)*0.96,max(means)*1.03); plt.tight_layout()
plt.savefig(os.path.join(OUT,'fig_wavelet_ablation.png'),dpi=150); plt.close()

# Fig: residual demand vs cuts (from committed CSV + cuts)
rd=os.path.join(ROOT,'results','residual_demand_monthly.csv')
import csv as _csv
mcs=[float(r['merit_cost']) for r in _csv.DictReader(open(rd))]
cuts=[3.0,5.8,81.7,81.7,5056.6,10000.0,10000.0,97.5,81.7,8.0,6.1,3.2,3.0,4.4,81.7,81.7,81.7,81.7]
fig,ax=plt.subplots(figsize=(7.8,4.2))
ax.step(range(1,13),mcs,where='mid',color=GREEN,lw=2.2,label='merit cost of monthly zone demand')
ax.plot(np.linspace(1,12,len(cuts)),cuts,'o--',color=NAVY,lw=1.4,ms=5,label='Plan4RES per-stage water values')
ax.set_yscale('log'); ax.set_xlabel('month / stage (rescaled)'); ax.set_ylabel('marginal value [$/MWh]')
ax.set_title('Stage trajectory: same merit-order levels; phasing differs (r=0.06)')
ax.legend(fontsize=8.5); plt.tight_layout()
plt.savefig(os.path.join(OUT,'fig_residual_demand.png'),dpi=150); plt.close()
print('figures written to', OUT, ':', sorted(os.listdir(OUT)))
