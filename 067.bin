from pathlib import Path
import re,json
ROOT=Path(__file__).resolve().parents[1]
s='GPAAGSSGTLELPSTLSLYKSGALDEAAAYQSRDYYNFPLALAGPPPPPPPPHPHARIKLENPLDYGSAWAAAAAQCRYGDLASLHGAGAAGPGSGSPSAAASSSWHTLFTAEEGQLYGP'
t=(ROOT/'background/bmr53105_3.str').read_text();ref=''.join(re.search(r'_Entity.Polymer_seq_one_letter_code\s*\n;\n(.*?)\n;',t,re.S).group(1).split())
assert s=='GP'+ref[:-1]
(ROOT/'background/WT_user_construct.fasta').write_text('>WT_Tau5_user_construct_120aa_expected_MW_12.1kDa\n'+s+'\n')
(ROOT/'background/BMRB53105_Tau5.fasta').write_text('>BMRB53105_Tau5_119aa\n'+ref+'\n')
(ROOT/'background/sequence_comparison.json').write_text(json.dumps(dict(user_length=len(s),bmrb_length=len(ref),relationship='user = GP + BMRB53105 sequence without terminal C',expected_monomer_mass_kDa=12.1,sequence_1_based_W_positions=[70,106],AR_W_positions=[397,433],mapping='User positions 3–120 match BMRB positions 1–118 and author AR positions 330–447 of BMRB51479.',source='https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr53105/bmr53105_3.str'),indent=2))
print('Validated exact sequence relationship:',len(s),len(ref))
