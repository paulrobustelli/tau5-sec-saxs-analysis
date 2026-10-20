from fit_experiments import *
for cutoff in [20,22,24]:
 ii=np.flatnonzero(r.Rg_CA_A>=cutoff);run(f'REST_WT_remove_Rg_lt{cutoff}',rr[:,ii],{'base_pool':'REST_WT','experiment':f'Remove Rg<{cutoff}A','source_indices':ii.tolist()},.1*len(r))
low=[];high=[];bins=np.floor(r.Rg_CA_A/2).astype(int)
for b in sorted(set(bins)):
 ii=np.flatnonzero(bins==b);ii=ii[np.argsort(r.loc[ii,'aromatic_hydrophobic_pairs_sep6_50'].to_numpy(),kind='stable')];half=len(ii)//2
 if half>=3:low.extend(ii[:half]);high.extend(ii[-half:])
for lab,ii in [('low_aromatic_Rgmatched',low),('high_aromatic_Rgmatched',high)]:run('REST_WT_'+lab,rr[:,ii],{'base_pool':'REST_WT','experiment':lab,'source_indices':list(map(int,ii))},.1*len(r))
for seed in [1,2]:
 rng=np.random.default_rng(seed);ii=np.sort(rng.choice(np.flatnonzero(r.Rg_CA_A<19),nadd,replace=False))
 for name in ['IDPForge_WT','bAIes_WT','ASTEROIDS_WT']:
  d,w,raw,qq,yy,ss=load(name);run(f'{name}_add_compact_REST_seed{seed}',np.c_[raw,rr[:,ii]],{'base_pool':name,'experiment':'add_compact_REST_replicate','seed':seed,'added_N':len(ii),'added_REST_indices':ii.tolist()},.1*len(d),donor_start=len(d))
