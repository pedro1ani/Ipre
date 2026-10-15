"""Construye matrices en mM y ratios; no entrena ni imputa. Python estándar."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import sqlite3

ROOT = Path(__file__).resolve().parent

def require(ok, message):
    if not ok:
        raise ValueError(message)

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def write_csv(path, rows, fields):
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'normalizacion')
    args = parser.parse_args()
    config_path = ROOT/'config_normalizacion.json'
    cfg = json.loads(config_path.read_text(encoding='utf-8'))
    inputs = [ROOT/'data/concentraciones_mM.db', ROOT/'sources/ecoli_auditada_v3.db', config_path, Path(__file__).resolve()]
    before = {str(p.relative_to(ROOT)): sha(p) for p in inputs}
    db = sqlite3.connect(inputs[0].resolve().as_uri()+'?mode=ro', uri=True)
    db.row_factory = sqlite3.Row
    source = sqlite3.connect(inputs[1].resolve().as_uri()+'?mode=ro', uri=True)
    source.row_factory = sqlite3.Row
    metabolites = list(cfg['ishii_analyte_ids'])
    targets = cfg['targets']
    records = [dict(r) for r in db.execute('SELECT * FROM metabolite_measurements_mM')]
    ishii = [r for r in records if r['dataset_id']=='ISHII2007' and r['compartment']=='intracellular' and r['metabolite_id'] in metabolites and r['analyte_id']==cfg['ishii_analyte_ids'][r['metabolite_id']]]
    metadata = {r['sample_id']: dict(r) for r in source.execute('SELECT * FROM ishii_samples')}
    available = sorted({r['sample_key'] for r in ishii})
    cohort = [s for s in available if metadata.get(s, {}).get('dilution_rate_h_inv')==0.2 or s in cfg['source_only_samples']]
    require(len(cohort)==30, 'Se esperaban 30 muestras Ishii')
    require('KO05x' in cohort and 'KO05' in cohort, 'Mantener ambos cultivos pfkA')
    references = {}
    for met in metabolites:
        selected = [r for r in records if r['dataset_id']=='MMC2' and r['condition_id']==cfg['reference_condition'] and r['metabolite_id']==met and r['compartment']=='intracellular']
        require(len(selected)==1, 'Referencia duplicada o ausente: '+met)
        r=selected[0]
        require(r['value_mM'] is not None and math.isfinite(r['value_mM']) and r['value_mM']>0, 'Referencia no positiva: '+met)
        require(r['absence_code']=='reported' and not r['is_imputed'] and not r['composite'], 'Referencia no admisible: '+met)
        require(r['uncertainty_type']=='SD' and r['uncertainty_mM'] is not None and r['uncertainty_mM']>=0 and r['source_uncertainty_cell'], 'SD sin documentar: '+met)
        references[met]=r
    conditions = {r['condition_id']:dict(r) for r in source.execute("SELECT * FROM conditions WHERE dataset_id='MMC2' AND condition_type='carbon_source'")}
    selected_mmc2=[]
    condition_rows=[]
    for condition in cfg['mmc2_conditions']:
        panel=[r for r in records if r['dataset_id']=='MMC2' and r['condition_id']==condition and r['metabolite_id'] in metabolites]
        require(condition in conditions and len(panel)==7 and len({r['metabolite_id'] for r in panel})==7, 'Panel MMC2 incompleto/duplicado')
        require(all(r['compartment']=='intracellular' and not r['composite'] and not r['is_imputed'] and r['mapping_status'] in ('exact_base_unique','confirmed_intracellular_source') for r in panel), 'Identidad MMC2 no admisible')
        selected_mmc2.extend(panel)
        condition_rows.append(dict(conditions[condition], admitted=True, observed=sum(r['value_mM'] is not None for r in panel), missing=sum(r['value_mM'] is None for r in panel), reference_row=condition==cfg['reference_condition'], reason='Identidad resuelta, intracelular, estacionario; ausencias se conservan. No implica equivalencia biológica con Ishii.', evaluation_role='auxiliary_training_only_not_validation_target'))
    selected=[r for r in ishii if r['sample_key'] in cohort]+selected_mmc2
    cells={}
    long=[]
    for r in selected:
        sample=r['sample_key'] if r['dataset_id']=='ISHII2007' else r['condition_id']
        met=r['metabolite_id'];ref=references[met];key=(sample,met)
        require(key not in cells, 'Celda duplicada: '+str(key))
        value=r['value_mM'];ratio=None if value is None else value/ref['value_mM']
        require(value is None or (math.isfinite(value) and value>=0 and math.isclose(ratio*ref['value_mM'],value,rel_tol=1e-12,abs_tol=1e-14)), 'Reversibilidad/valor inválido')
        cells[key]=(value,ratio)
        long.append(dict(sample_id=sample,dataset_id=r['dataset_id'],metabolite_id=met,role='target' if met in targets else 'predictor',value_mM=value,ratio=ratio,observed=value is not None,absence_code=r['absence_code'],measurement_uid=r['measurement_uid'],analyte_id=r['analyte_id'],protocol=r['protocol'],source_file=r['source_file'],source_sheet=r['source_sheet'],source_cell=r['source_cell'],source_sha256=r['source_file_sha256'],uncertainty_mM=r['uncertainty_mM'],uncertainty_type=r['uncertainty_type'],reference_uid=ref['measurement_uid'],reference_mM=ref['value_mM'],reference_SD_mM=ref['uncertainty_mM'],is_reference_cell=r['measurement_uid']==ref['measurement_uid'],ratio_uncertainty_status='not_propagated_shared_denominator',is_imputed=False))
    require(len(cells)==266,'Se esperaban 266 celdas')
    require(all(cells[(cfg['reference_condition'],m)][1]==1 for m in metabolites),'Referencia normalizada distinta de 1')
    sample_rows=[]
    for s in available:
        md=metadata.get(s,{})
        sample_rows.append(dict(sample_id=s,dataset_id='ISHII2007',included=s in cohort,perturbation_group='pfkA' if s in ('KO05','KO05x') else (md.get('knockout_gene') or 'WT'),dilution_h_inv=md.get('dilution_rate_h_inv'),carbon_source='glucose',protocol_family='mixed_analyte_protocols',metadata_status='source_only_culture_no_S2_row' if s=='KO05x' else 'documented',reason='retain_previous_cohort_pending_metadata' if s=='KO05x' else ('dilution_0.2' if s in cohort else 'other_dilution_excluded')))
    for r in condition_rows:
        sample_rows.append(dict(sample_id=r['condition_id'],dataset_id='MMC2',included=True,perturbation_group=r['condition_id'],dilution_h_inv=None,carbon_source=r['carbon_source'],protocol_family='Gerosa',metadata_status='documented',reason='auxiliary_training_condition'))
    summary={}
    out=args.output.resolve()
    out.mkdir(parents=True,exist_ok=False)
    for scenario,samples in [('ishii',cohort),('ishii_mmc2',cohort+cfg['mmc2_conditions'])]:
        for kind,idx in [('mM',0),('R',1),('observado',None)]:
            matrix=[dict(sample_id=s,**{m:(int(cells[(s,m)][0] is not None) if idx is None else cells[(s,m)][idx]) for m in metabolites}) for s in samples]
            write_csv(out/f'matriz_{scenario}_{kind}.csv',matrix,['sample_id']+metabolites)
        observed=sum(cells[(s,m)][0] is not None for s in samples for m in metabolites)
        target_observed=sum(cells[(s,m)][0] is not None for s in samples for m in targets)
        summary[scenario]=dict(rows=len(samples),columns=7,observed=observed,missing=len(samples)*7-observed,target_observed=target_observed,target_missing=len(samples)*4-target_observed)
    write_csv(out/'trazabilidad_celdas.csv',long,list(long[0]))
    write_csv(out/'muestras_y_exclusiones.csv',sample_rows,list(sample_rows[0]))
    write_csv(out/'condiciones_mmc2.csv',condition_rows,list(condition_rows[0]))
    ref_rows=[dict(metabolite_id=m,reference_condition=cfg['reference_condition'],value_mM=r['value_mM'],SD_mM=r['uncertainty_mM'],uncertainty_type=r['uncertainty_type'],measurement_uid=r['measurement_uid'],source_file=r['source_file'],source_sheet=r['source_sheet'],source_cell=r['source_cell'],source_SD_cell=r['source_uncertainty_cell'],source_sha256=r['source_file_sha256']) for m,r in references.items()]
    write_csv(out/'referencias_glucosa.csv',ref_rows,list(ref_rows[0]))
    require(before=={str(p.relative_to(ROOT)):sha(p) for p in inputs},'Entradas modificadas')
    manifest=dict(status='SUCCESS',stage='normalized_no_imputation',summary=summary,input_sha256=before,checks=dict(positive_observed_references=True,reference_SD_documented=True,unique_cells=True,missing_preserved=True,reversible=True,glucose_ratios_one=True,inputs_unchanged=True),output_sha256={p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()})
    (out/'manifest_normalizacion.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps(summary,indent=2))
    db.close();source.close()

if __name__=='__main__':
    main()
