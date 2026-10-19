"""Gráficos de resultados existentes; no vuelve a entrenar ni imputar."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from imputadores import METHODS
ROOT = Path(__file__).resolve().parent

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--results', type=Path, default=ROOT/'resultados/ishii')
    out = parser.parse_args().results
    pred = pd.read_csv(out/'validacion/predicciones.csv')
    figures = out/'figuras'
    figures.mkdir(exist_ok=True)
    # Una figura por prueba: medición original ocultada frente a estimación.
    for experiment, part in pred.groupby('experiment'):
        metabolites = list(part.metabolite_id.unique())
        fig, axes = plt.subplots(len(metabolites), 3, figsize=(12, 3*len(metabolites)), squeeze=False)
        for i, met in enumerate(metabolites):
            for j, method in enumerate(METHODS):
                ax = axes[i,j]
                values = part[(part.metabolite_id==met)&(part.method==method)].dropna(subset=['estimate_mM'])
                supported = values.support=='supported'
                for selected, color, label in [(supported, '#147d92', 'Con soporte'), (~supported, '#d78b35', 'General o limitado')]:
                    sub = values[selected]
                    ax.scatter(sub.truth_mM, sub.estimate_mM, color=color, s=23, alpha=.7, label=label)
                maximum = max(values.truth_mM.max(), values.estimate_mM.max(), .001) if len(values) else 1
                ax.plot([0,maximum], [0,maximum], '--', color='gray', lw=1)
                ax.set_title(f'{method}\n{met}; n={len(values)}', fontsize=10)
                ax.set_xlabel('Medición ocultada (mM)')
                ax.set_ylabel('Estimación (mM)')
                ax.grid(alpha=.15)
        axes[0,0].legend(fontsize=8)
        fig.suptitle(f'{experiment}: valor conocido vs imputado\nTodos los candidatos numéricos; diagonal = predicción exacta')
        fig.tight_layout(rect=(0,0,1,.95))
        fig.savefig(figures/f'{experiment}_real_vs_imputado.png', dpi=150)
        plt.close(fig)
    print(f'Gráficos: {figures}')

if __name__=='__main__':
    main()
