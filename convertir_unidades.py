"""Reconstruye una capa de concentraciones en mM sin alterar fuentes ni imputar."""
from pathlib import Path
import json,hashlib,sqlite3,sys
import pandas as pd
import numpy as np
import openpyxl
from pypdf import PdfReader
R=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    out=R/'data';out.mkdir(exist_ok=True)
    db=out/'concentraciones_mM.db'
    if db.exists():raise FileExistsError('Use una copia con data vacía para reconstruir; no se sobrescribe SQLite.')
    before={p.name:sha(p) for p in (R/'sources').iterdir() if p.is_file()}
    w2=openpyxl.load_workbook(R/'sources/mmc2.xlsx',data_only=True)
    w3=openpyxl.load_workbook(R/'sources/mmc3.xlsx',data_only=True)
    wi=openpyxl.load_workbook(R/'sources/ishii.xlsx',data_only=True)
    cf=w2['Conversion Factors'];V=float(cf['B6'].value)
    assert V==2.3 and cf['B5'].value=='All conditions' and 'ml/gCDW' in cf['A6'].value
    methods=PdfReader(R/'sources/Gerosa_methods.pdf').pages[8].extract_text()
    assert 'intracellular cytosolic volume' in methods and 'µmol/gCDW' in methods
    factor=dict(factor_id='GEROSA_CYTOSOLIC_VOLUME',value=V,unit='mL/gCDW',definition='intracellular cytosolic volume per dry biomass',scope='All conditions; MMC2 and MMC3 intracellular metabolites',source_file='mmc2.xlsx',source_sheet='Conversion Factors',source_cell='B6',definition_cells='A5:B6',methods_source='Gerosa_methods.pdf',methods_pdf_page=9,methods_section='Intracellular metabolite concentrations by LC-MS/MS',factor_sd=None,factor_uncertainty_status='not_reported',origin_limit='Source says previously determined; original measurement protocol and uncertainty not specified here')
    source=R/'sources/ecoli_auditada_v3.db'
    with sqlite3.connect(source.resolve().as_uri()+'?mode=ro',uri=True) as c:
        d=pd.read_sql_query("SELECT * FROM audit_measurements_readable WHERE entity_type='metabolite'",c)
    assert len(d)==20483 and not d.measurement_uid.duplicated().any()
    d['source_file']=d.dataset_id.map({'MMC2':'mmc2.xlsx','MMC3':'mmc3.xlsx','ISHII2007':'ishii.xlsx'})
    d['source_file_sha256']=d.source_file.map(before)
    # Verificar cada valor no nulo contra su celda de origen, no contra resultados imputados.
    books={'MMC2':w2,'MMC3':w3,'ISHII2007':wi}
    for row in d.itertuples():
        raw=books[row.dataset_id][row.source_sheet][row.source_cell].value
        if pd.notna(row.value_original):assert isinstance(raw,(int,float)) and np.isclose(raw,row.value_original,rtol=1e-12,atol=1e-14),(row.measurement_uid,raw)
        if row.dataset_id in ['MMC2','MMC3'] and pd.notna(row.uncertainty_original):
            sd=books[row.dataset_id][row.source_sheet][row.source_uncertainty_cell].value
            assert np.isclose(sd,row.uncertainty_original,rtol=1e-12,atol=1e-14)
    ger=d.dataset_id.isin(['MMC2','MMC3'])
    assert set(d.loc[ger,'unit_original'])=={'umol gCDW^-1'}
    assert set(d.loc[~ger,'unit_original'])=={'mM'}
    d['conversion_divisor']=np.where(ger,V,1.)
    d['factor_id']=np.where(ger,factor['factor_id'],'IDENTITY_mM')
    d['intracellular_volume_mL_gCDW']=np.where(ger,V,np.nan)
    d['unit_final']='mM'
    d['value_mM']=d.value_original/d.conversion_divisor
    d['uncertainty_mM']=d.uncertainty_original/d.conversion_divisor
    d['conversion_formula']=np.where(ger,'value_original / 2.3','value_original (unchanged)')
    d['conversion_uncertainty_note']=np.where(ger,'SD scaled with fixed V; excludes unreported uncertainty of V','Original uncertainty retained without adding estimates')
    d['processing_stage']='unit_conversion_only_no_reference_normalization_no_imputation'
    assert d.value_mM.isna().equals(d.value_original.isna())
    assert np.allclose(d.value_mM*d.conversion_divisor,d.value_original,equal_nan=True)
    assert np.array_equal(d.loc[~ger,'value_mM'].fillna(-999),d.loc[~ger,'value_original'].fillna(-999))
    assert (d.is_imputed==0).all()
    with sqlite3.connect(db) as c:
        d.to_sql('metabolite_measurements_mM',c,index=False)
        pd.DataFrame([factor]).to_sql('conversion_factors',c,index=False)
        c.executescript("CREATE UNIQUE INDEX measurement_uid ON metabolite_measurements_mM(measurement_uid); CREATE VIEW intracellular_mM AS SELECT * FROM metabolite_measurements_mM WHERE compartment='intracellular'; CREATE VIEW extracellular_mM AS SELECT * FROM metabolite_measurements_mM WHERE compartment='extracellular';")
        assert c.execute('pragma integrity_check').fetchone()[0]=='ok'
    for key,g in d.groupby('dataset_id'):
        g.to_csv(out/(key+'_mM.csv'),index=False,encoding='utf-8-sig',float_format='%.17g')
    pd.DataFrame([factor]).to_csv(out/'conversion_factors.csv',index=False,encoding='utf-8-sig')
    (out/'factor_evidence.txt').write_text('MMC2 Conversion Factors A5:B6: '+str([[cf.cell(i,j).value for j in [1,2]] for i in [5,6]])+'\n\nGerosa methods, PDF page 9:\n'+methods)
    summary=d.groupby(['dataset_id','compartment'],as_index=False).agg(cells=('measurement_uid','size'),observed=('value_mM','count'))
    summary['missing']=summary.cells-summary.observed;summary.to_csv(out/'summary.csv',index=False)
    assert before=={p.name:sha(p) for p in (R/'sources').iterdir() if p.is_file()}
    manifest=dict(source_sha256=before,script_sha256=sha(Path(__file__)),factor=factor,checks={'all_original_values_checked_against_workbooks':True,'Gerosa_SD_checked':True,'missing_pattern_preserved':True,'Ishii_unchanged':True,'sources_unchanged':True,'reverse_conversion_passed':True},python=sys.version,output_sha256={p.name:sha(p) for p in out.iterdir() if p.is_file()})
    (R/'manifest.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False))
    print(summary.to_string(index=False))
if __name__=='__main__':main()
