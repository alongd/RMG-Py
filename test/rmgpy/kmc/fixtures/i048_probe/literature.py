"""Arithmetic cross-check of explicitly identified original conformation inputs.

Command: rmg_env python .../i048_probe/literature.py
These are publication inputs from indexed primary excerpts, not new QM data.
"""
import math
from common import SCRATCH,save

R=8.314472
DATA={
    'retrieval_limit':'Complete Yoon 1975 and Williams/Flory 1969 PDFs could not be fetched. Original-paper indexed excerpts and publisher/university metadata were reached. Full RIS matrix reconstruction and an absolute entropy-per-dyad comparison are not asserted.',
    'Yoon_1975':{
        'citation':'D. Y. Yoon, P. R. Sundararajan and P. J. Flory, Conformational Characteristics of Polystyrene, Macromolecules 8 (1975) 776–783',
        'doi':'https://doi.org/10.1021/ma60048a019',
        'indexed_primary_pdf':'https://electronicsandbooks.com/edt/manual/Magazine/M/Macromolecules/1975%20%28Vol%208%29/No06%28691-959%29/776.pdf',
        'location':'p. 781, rounded prefactors and equations (7)–(9)',
        'eta_prefactor':0.8,'omega_prefactor':1.3,
        'R_ln_eta_prefactor_J_mol_K':R*math.log(.8),
        'R_ln_omega_prefactor_J_mol_K':R*math.log(1.3),
        'interpretation':'relative statistical-weight prefactors give local well-shape entropy differences; these are not a complete entropy increment or a QM-minus-GAV correction'},
    'Williams_Flory_1969':{
        'citation':'A. D. Williams and P. J. Flory, Stereochemical equilibrium and configurational statistics in polystyrene and its oligomers, J. Am. Chem. Soc. 91 (1969) 3111–3118',
        'doi':'https://doi.org/10.1021/ja01040a001',
        'indexed_primary_pdf':'https://www.electronicsandbooks.com/edt/manual/Magazine/J/Journal%20of%20the%20American%20Chemical%20Society%20US/1969%20%20%28vol%20091%29/12%20%20%283111-3408%29/3111-3118.pdf',
        'location':'discussion following equation (28)',
        'relative_H_cal_mol':-700.,'relative_S_cal_mol_K':-1.4,
        'relative_H_kJ_mol':-.700*4.184,'relative_S_J_mol_K':-1.4*4.184,
        'interpretation':'reported relative conformer preference with entropy opposing stabilization; this is not the atactic propagation correction'},
    'Khare_Paulaitis_1994':{
        'citation':'R. Khare and M. E. Paulaitis, Molecular simulations of cooperative ring flip motions in single chains of polystyrene, Chemical Engineering Science 49 (1994) 2867–2879',
        'doi':'https://doi.org/10.1016/0009-2509(94)E0105-Y',
        'accessible_primary_abstract':'https://pure.johnshopkins.edu/en/publications/molecular-simulations-of-cooperative-ring-flip-motions-in-single-',
        'interpretation':'hexamer calculations explicitly couple phenyl and backbone torsions; the reached abstract provides no absolute entropy-per-dyad benchmark'},
    'fit_to_literature':False,
}

if __name__=='__main__':
    save(SCRATCH/'literature_crosscheck.json',DATA)
    print('I048 literature arithmetic: eta prefactor %.6f J/mol/K; omega prefactor %.6f J/mol/K; reported conformer difference %.6f kJ/mol and %.6f J/mol/K'%(
        DATA['Yoon_1975']['R_ln_eta_prefactor_J_mol_K'],DATA['Yoon_1975']['R_ln_omega_prefactor_J_mol_K'],
        DATA['Williams_Flory_1969']['relative_H_kJ_mol'],DATA['Williams_Flory_1969']['relative_S_J_mol_K']))
    print('I048 literature access limit: '+DATA['retrieval_limit'])
