"""Descriptivos y gráficos exploratorios de valores observados; sin imputación."""
from pathlib import Path
import argparse,json,hashlib,importlib.metadata
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
ROOT=Path(__file__).resolve().parent
LABELS=['Fumarato','PEP','s7p','Succinato','Piruvato','α-cetoglutarato','Citrato']

def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'descriptivos');a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
 paths=[ROOT/'normalizacion'/f'matriz_ishii_mmc2_{unit}.csv' for unit in ['mM','R']]
 frames={unit:pd.read_csv(path,index_col=0) for unit,path in zip(['mM','R'],paths)}
 groups={'Ishii':lambda idx:~idx.str.startswith('MMC2_'),'MMC2':lambda idx:idx.str.startswith('MMC2_'),'Combinado':lambda idx:np.ones(len(idx),dtype=bool)}
 stats=[]
 for unit,df in frames.items():
  for group,select in groups.items():
   d=df.loc[select(df.index)]
   for met in d:
    x=d[met].dropna();n=len(x);mean=x.mean();sd=x.std(ddof=1);q1=x.quantile(.25);q3=x.quantile(.75)
    stats.append(dict(group=group,unit=unit,metabolite=met,total=len(d),n=n,missing=len(d)-n,coverage_pct=100*n/len(d),mean=mean,median=x.median(),SD=sd,Q1=q1,Q3=q3,IQR=q3-q1,CV_pct=100*sd/mean if mean!=0 else np.nan,min=x.min(),max=x.max()))
 pd.DataFrame(stats).to_csv(out/'estadisticos.csv',index=False,encoding='utf-8-sig')
 ratios=frames['R'];ish=ratios.loc[groups['Ishii'](ratios.index)];ger=ratios.loc[groups['MMC2'](ratios.index)]
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
 figures=[]
 fig,axs=plt.subplots(2,4,figsize=(14,8));rng=np.random.default_rng(20261015)
 for j,(met,label) in enumerate(zip(ratios.columns,LABELS)):
  ax=axs.flat[j]
  for pos,(name,data,color) in enumerate([('Ishii',ish,'#196b9a'),('MMC2',ger,'#d47929')],1):
   x=data[met].dropna();ax.boxplot(x,positions=[pos],widths=.45,showfliers=False,patch_artist=True,boxprops={'facecolor':color,'alpha':.22},medianprops={'color':'black'})
   ax.scatter(pos+rng.uniform(-.13,.13,len(x)),x,s=25,color=color,alpha=.75)
  ax.axhline(1,color='#999999',ls='--',lw=1);ax.set_yscale('log');ax.set_xticks([1,2],['Ishii\nn='+str(ish[met].count()),'MMC2\nn='+str(ger[met].count())]);ax.set_title(label);ax.set_ylabel('Ratio respecto a glucosa');ax.grid(axis='y',alpha=.15)
 axs.flat[7].axis('off');axs.flat[7].text(0,.85,'Solo valores observados.\nEje vertical logarítmico.\nLínea punteada: referencia = 1.\nPuntos: muestras/condiciones.\nCajas: Q1–Q3 y mediana.\nBigotes: hasta 1,5 IQR.\nEscalas distintas entre metabolitos.\nMMC2 incluye la fila de referencia.',va='top',linespacing=1.6)
 fig.suptitle('Distribuciones antes de imputar',fontsize=18);fig.tight_layout(rect=(0,0,1,.95));figures.append(('01_distribuciones',fig))
 fig,ax=plt.subplots(figsize=(12,11));vals=ratios.to_numpy();masked=np.ma.masked_invalid(np.log2(np.where(vals>0,vals,np.nan)));cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#d6d6d6');lim=max(abs(masked.min()),abs(masked.max()));im=ax.imshow(masked,aspect='auto',cmap=cmap,vmin=-lim,vmax=lim)
 ax.set_xticks(range(7),LABELS,rotation=25,ha='right');ax.set_yticks(range(len(ratios)),ratios.index,fontsize=8);ax.axhline(29.5,color='black',lw=2);fig.colorbar(im,ax=ax,label='log₂(ratio): 0 = referencia; +1 = doble; −1 = mitad');ax.set_title('Patrones por muestra y metabolito\nGris = faltante; línea negra separa Ishii y MMC2',fontsize=16);fig.tight_layout();figures.append(('02_mapa_muestras',fig))
 fig,axs=plt.subplots(1,2,figsize=(13,5));x=np.arange(7);width=.35
 for k,(name,df,color) in enumerate([('Ishii',ish,'#196b9a'),('MMC2',ger,'#d47929')]):
  cov=df.notna().mean()*100;bars=axs[0].bar(x+(k-.5)*width,cov,width,label=name,color=color)
  for b,n in zip(bars,df.count()):axs[0].text(b.get_x()+b.get_width()/2,b.get_height()+2,str(n),ha='center',fontsize=8)
 axs[0].set_xticks(x,LABELS,rotation=35,ha='right');axs[0].set_ylim(0,115);axs[0].set_ylabel('Cobertura (%)');axs[0].legend();axs[0].set_title('Cobertura por metabolito\nEtiquetas = número de valores medidos')
 count=ish.notna().sum(axis=1);freq=count.value_counts().reindex(range(8),fill_value=0);axs[1].bar(freq.index,freq,color='#196b9a');axs[1].set_xticks(range(8));axs[1].set_xlabel('Metabolitos observados de 7');axs[1].set_ylabel('Número de muestras Ishii');axs[1].set_title('Información disponible por muestra');fig.tight_layout();figures.append(('03_cobertura',fig))
 # Correlación por pares en Ishii, sin rellenar faltantes; mínimo 8 pares.
 n=ish.notna().astype(int).T.dot(ish.notna().astype(int));corr=ish.corr(method='spearman',min_periods=8);corr.to_csv(out/'correlacion_spearman_ishii.csv',encoding='utf-8-sig');n.to_csv(out/'n_pares_ishii.csv',encoding='utf-8-sig')
 fig,ax=plt.subplots(figsize=(10,8));im=ax.imshow(np.ma.masked_invalid(corr),vmin=-1,vmax=1,cmap='RdBu_r');ax.set_xticks(range(7),LABELS,rotation=35,ha='right');ax.set_yticks(range(7),LABELS)
 for i in range(7):
  for j in range(7):
   value=corr.iloc[i,j];ax.text(j,i,('—' if pd.isna(value) else f'{value:.2f}')+f'\nn={n.iloc[i,j]}',ha='center',va='center',color='white' if pd.notna(value) and abs(value)>.65 else 'black',fontsize=10)
 fig.colorbar(im,ax=ax,label='ρ de Spearman');ax.set_title('Asociaciones exploratorias dentro de Ishii\nPares observados; mínimo 8 pares; sin imputación',fontsize=15);fig.tight_layout();figures.append(('04_correlaciones_ishii',fig))
 with PdfPages(out/'graficos_exploratorios.pdf') as pdf:
  for name,fig in figures:fig.savefig(out/(name+'.png'),dpi=180);pdf.savefig(fig);plt.close(fig)
 count.rename('observed_metabolites').to_csv(out/'cobertura_muestras_ishii.csv',encoding='utf-8-sig')
 manifest={'stage':'exploratory_observed_only','input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'versions':{m:importlib.metadata.version(m) for m in ['numpy','pandas','matplotlib','scipy']},'SD_ddof':1,'quantiles':'linear','correlation':'Spearman pairwise complete, min 8','jitter_seed':20261015,'outputs_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir()}}
 (out/'manifest_descriptivos.json').write_text(json.dumps(manifest,indent=2)+'\n')
 print(pd.DataFrame(stats).query("unit=='R' and group!='Combinado'")[['group','metabolite','n','median','CV_pct','min','max']].to_string(index=False));print('\nCORRELACIONES\n',corr.round(3).to_string());print('\nMUESTRAS CON <=2 VALORES\n',count[count<=2].to_string())
if __name__=='__main__':main()
