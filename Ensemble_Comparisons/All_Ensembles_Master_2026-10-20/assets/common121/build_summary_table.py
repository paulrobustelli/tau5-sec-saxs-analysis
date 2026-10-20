"""One table of every available stage, with missing results explicitly absent."""
from pathlib import Path
import pandas as pd
O=Path(__file__).resolve().parent
s=pd.read_csv(O/'structural_comparison.csv');e=pd.read_csv(O/'experimental_agreement.csv');e['column']=e.apply(lambda r:'SAXS_chi2_mean' if r.datatype=='SAXS' else 'CS_'+r.datatype+'_RMSE_ppm',axis=1)
v=e.pivot(index=['ensemble','stage'],columns='column',values='value').reset_index();n=e.pivot(index=['ensemble','stage'],columns='datatype',values='N_observations').add_prefix('N_observations_').reset_index();table=s.merge(v,on=['ensemble','stage'],how='left',validate='one_to_one').merge(n,on=['ensemble','stage'],how='left',validate='one_to_one');table['normalized_Kish']=table.Neff/table.N;table.to_csv(O/'master_comparison_table.csv',index=False)
print(table[['ensemble','stage','Rg_CA_A','R2_main_H_percent','R2_minor_H_percent','R3_H_percent','SAXS_chi2_mean']].to_string(index=False))
