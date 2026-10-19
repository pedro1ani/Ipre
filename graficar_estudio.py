"""Figuras de validación y candidatos; no cambia modelos ni estimaciones."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
ROOT=Path(__file__).resolve().parent
METHODS=['SimpleImputer','KNNImputer','IterativeImputer'];LABELS=['Mediana','KNN','Iterativo'];SC=['ishii','ishii_mmc2'];T=['fum[c]','pep[c]','s7p[c]','succ[c]'];TL=['Fumarato','PEP','s7p','Succinato']

def main():
 p=argparse.ArgumentParser();p.add_argument('--results',type=Path,default=ROOT/'resultados_estudio');a=p.parse_args();out=a.results/'figuras';out.mkdir(exist_ok=True)
 metrics=pd.read_csv(a.results/'validacion/metricas.csv');plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
 figs=[];fig,axs=plt.subplots(1,4,figsize=(15,5))
 for j,(met,label) in enumerate(zip(T,TL)):
  ax=axs[j]
  for k,(sc,color) in enumerate(zip(SC,['#196b9a','#d47929'])):
   q=metrics[(metrics.policy=='available')&(metrics.scope=='common_all_six')&(metrics.experiment=='perturbation')&(metrics.scenario==sc)&(metrics.metabolite_id==met)].set_index('method').loc[METHODS]
   ax.bar(np.arange(3)+(k-.5)*.35,q.MAE_mM,.35,color=color,label='Ishii' if k==0 else 'Ishii + MMC2')
  ax.set_xticks(range(3),LABELS,rotation=25,ha='right');ax.set_title(label+f' · n={int(q.predictions.iloc[0])}');ax.set_ylabel('MAE (mM)');ax.grid(axis='y',alpha=.15)
 axs[0].legend(fontsize=9);fig.suptitle('Exclusión por perturbación · mismas celdas en los seis casos\nTodos los candidatos numéricos, incluidos fallbacks y ajustes sin convergencia',fontsize=14);fig.tight_layout(rect=(0,0,1,.87));figs.append(('01_MAE_perturbacion',fig))
 sub=metrics[(metrics.policy=='supported_plus_baseline')&(metrics.scope=='own')];experiments=['perturbation','real_patterns','analytical_series'];expnames=['Perturbación','Patrones','Serie']
 values=[];texts=[]
 for sc in SC:
  for method in METHODS:
   rr=[];tt=[]
   for exp in experiments:
    for met in T:
     row=sub[(sub.scenario==sc)&(sub.method==method)&(sub.experiment==exp)&(sub.metabolite_id==met)].iloc[0];rr.append(row.coverage_pct);tt.append(f'{int(row.predictions)}/{int(row.attempts)}')
   values.append(rr);texts.append(tt)
 fig,ax=plt.subplots(figsize=(16,6));im=ax.imshow(values,vmin=0,vmax=100,cmap='YlGnBu',aspect='auto')
 for i in range(6):
  for j in range(12):ax.text(j,i,texts[i][j],ha='center',va='center',color='white' if values[i][j]>65 else 'black',fontsize=9)
 ax.set_xticks(range(12),[f'{e}\n{t}' for e in expnames for t in TL],rotation=30,ha='right');ax.set_yticks(range(6),[f'{s} · {m}' for s in ['Ishii','Ishii + MMC2'] for m in LABELS]);ax.axvline(3.5,color='white',lw=2);ax.axvline(7.5,color='white',lw=2);ax.axhline(2.5,color='white',lw=2);fig.colorbar(im,ax=ax,label='Cobertura con soporte + referencia mediana (%)');ax.set_title('Valores aceptados / intentos\nMediana se conserva como referencia general; KNN e iterativo requieren soporte',fontsize=14);fig.tight_layout();figs.append(('02_cobertura_operativa',fig))
 original=pd.read_csv(ROOT/'normalizacion/matriz_ishii_mM.csv',index_col=0);rng=np.random.default_rng(20261019)
 fig,axs=plt.subplots(1,2,figsize=(15,6))
 names=['Observados']+[f'{m}\n{s}' for s in ['Ishii','+ MMC2'] for m in LABELS]
 for ax,met,label in zip(axs,T[:2],TL[:2]):
  groups=[original[met]]
  for sc in SC:
   for method in METHODS:groups.append(pd.read_csv(a.results/'imputacion'/f'available_{sc}_{method}_mM.csv',index_col=0)[met])
  for k,series in enumerate(groups):
   finite=series.dropna();ax.boxplot(finite,positions=[k],widths=.48,showfliers=False,medianprops={'color':'black'})
   measured=original[met].notna().reindex(series.index)
   ax.scatter(k+rng.uniform(-.12,.12,int(measured.sum())),series[measured],color='#7895a5',s=15,alpha=.6)
   filled=series.notna()&~measured
   ax.scatter(k+rng.uniform(-.12,.12,int(filled.sum())),series[filled],color='#d47929',s=35,marker='D')
  ax.set_xticks(range(7),names,rotation=35,ha='right');ax.set_ylabel('Concentración (mM)');ax.set_title(label);ax.grid(axis='y',alpha=.15)
 fig.suptitle('Antes/después en las 30 muestras Ishii · todos los candidatos disponibles\nGris: medido; naranja: estimado. Incluye fallbacks y ajustes sin convergencia.',fontsize=14);fig.tight_layout(rect=(0,0,1,.87));figs.append(('03_antes_despues',fig))
 with PdfPages(out/'figuras_estudio.pdf') as pdf:
  for name,fig in figs:fig.savefig(out/(name+'.png'),dpi=180,bbox_inches='tight');pdf.savefig(fig,bbox_inches='tight');plt.close(fig)
 print('Tres figuras guardadas')
if __name__=='__main__':main()
