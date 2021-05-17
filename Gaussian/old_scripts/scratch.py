old_negs = {
    "mol_NN_cp_frag_cation1": 10,
    "mol_NN_cp_frag_cation2": 10,
    "mol_NN_cp_frag_neutral": 12,
    "mol_NN_ma_frag_neutral": 14,
    "mol_NN_pa_frag_cation1": 6,
    "mol_NN_pa_frag_cation2": 6,
    "mol_NN_pa_frag_neutral": 6,
    "mol_bp_NN_frag_cation1": 4,
    "mol_bp_NN_frag_cation2": 8,
    "mol_bp_NN_frag_neutral": 3,
    "mol_bp_cp_half_cation1": 29,
    "mol_bp_cp_half_neutral": 26,
    "mol_bp_cp_whol_cation2": 21,
    "mol_bp_cp_whol_neutral": 23,
    "mol_bp_ma_half_cation2": 17,
    "mol_bp_ma_half_neutral": 14,
    "mol_bp_ma_whol_cation2": 42,
    "mol_bp_ma_whol_neutral": 40,
    "mol_bp_pa_half_cation1": 9,
    "mol_bp_pa_half_cation2": 9,
    "mol_bp_pa_half_neutral": 9,
    "mol_bp_pa_whol_neutral": 60,
    "mol_ta_NN_frag_cation1": 10,
    "mol_ta_NN_frag_cation2": 9,
    "mol_ta_NN_frag_neutral": 6,
    "mol_ta_cp_half_cation2": 56,
    "mol_ta_cp_half_neutral": 56,
    "mol_ta_cp_whol_cation2": 20,
    "mol_ta_cp_whol_neutral": 22,
    "mol_ta_ma_half_cation2": 26,
    "mol_ta_ma_half_neutral": 24,
    "mol_ta_ma_whol_cation2": 18,
    "mol_ta_ma_whol_neutral": 19,
    "mol_ta_pa_half_cation1": 9,
    "mol_ta_pa_half_cation2": 9,
    "mol_ta_pa_half_neutral": 9
}

new_negs = {
  "mol_bp_cp_half_neutral": 26,
  "mol_bp_pa_half_neutral": 9,
  "mol_bp_pa_whol_neutral": 0,
  "mol_bp_cp_whol_neutral": 23,
  "mol_NN_ma_frag_neutral": 13,
  "mol_ta_ma_half_neutral": 26,
  "mol_NN_pa_frag_neutral": 6,
  "mol_ta_cp_half_neutral": 54,
  "mol_ta_pa_half_neutral": 9,
  "mol_ta_NN_frag_neutral": 6,
  "mol_bp_ma_half_neutral": 17,
  "mol_bp_NN_frag_neutral": 3,
  "mol_NN_cp_frag_neutral": 12
}

worse = []
better = []
same = []

for mol in new_negs.keys():
    if old_negs[mol] > new_negs[mol]:
        better.append(mol)
    elif old_negs[mol] < new_negs[mol]:
        worse.append(mol)
    elif old_negs[mol] == new_negs[mol]:
        same.append(mol)

print(worse)
print(better)
print(same)
