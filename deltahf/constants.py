"""Physical constants and model definitions shared across deltahf modules."""

HARTREE_TO_KCAL = 627.5094740631
HARTREE_TO_EV = 27.211386245988

# Atom type names for each classification scheme.
# All lists include F, S, Cl params; dynamic parameter filtering in atom_equivalents.py
# automatically drops any params with zero training examples (preserving CHNO-only behaviour).
PARAM_NAMES_4 = ["C", "H", "N", "O", "F", "S", "Cl"]
PARAM_NAMES_7 = ["C", "H", "N", "O", "C_prime", "N_prime", "O_prime", "F", "S", "S_prime", "Cl"]
PARAM_NAMES_HYBRID = [
    "C_sp3", "C_sp2", "C_sp", "H", "N_sp3", "N_sp2", "N_sp", "O_sp3", "O_sp2", "O_sp",
    "F_sp3", "S_sp3", "S_sp2", "S_sp", "Cl_sp3",
]
PARAM_NAMES_EXTENDED = [
    "C_sp3_3H", "C_sp3_2H", "C_sp3_1H", "C_sp3_0H",
    "C_sp2_2H", "C_sp2_1H", "C_sp2_0H", "C_sp",
    "H", "N_sp3", "N_sp2", "N_sp", "O_sp3", "O_sp2", "O_sp",
    "F", "S_sp3", "S_sp2", "Cl",
]
PARAM_NAMES_NEIGHBOUR = [
    "C_sp3_3H", "C_sp3_2H", "C_sp3_1H", "C_sp3_0H",
    "C_sp2_2H", "C_sp2_1H", "C_sp2_0H", "C_sp",
    "H",
    "N_sp3_C", "N_sp3_N", "N_sp3_O",
    "N_sp2_C", "N_sp2_N", "N_sp2_O",
    "N_sp_C",  "N_sp_N",  "N_sp_O",
    "O_sp3_C", "O_sp3_N", "O_sp3_O",
    "O_sp2_C", "O_sp2_N", "O_sp2_O",
    "O_sp_C",  "O_sp_N",  "O_sp_O",
    "F", "S_sp3", "S_sp2", "Cl",
]

PARAM_NAMES_BONDORDER = ["C_1", "C_2", "C_3", "H", "N_1", "N_2", "N_3", "O_1", "O_2", "O_3",
                          "F_1", "S_1", "S_2", "S_3", "Cl_1"]
PARAM_NAMES_BONDORDER_EXT = [
    "C_1_3H", "C_1_2H", "C_1_1H", "C_1_0H",
    "C_2_2H", "C_2_1H", "C_2_0H",
    "C_3_1H", "C_3_0H",
    "H",
    "N_1", "N_2", "N_3",
    "O_1", "O_2", "O_3",
    "F_1", "S_1", "S_2", "S_3", "Cl_1",
]
PARAM_NAMES_BONDORDER_AR = [
    "C_1", "C_ar", "C_2", "C_3", "H", "N_1", "N_ar", "N_2", "N_3", "O_1", "O_ar", "O_2", "O_3",
    "F_1", "S_1", "S_ar", "S_2", "S_3", "Cl_1",
]

# Maps model name -> (param_names, MoleculeResult atom_counts field)
MODEL_DEFS = {
    "4param":         (PARAM_NAMES_4,            "atom_counts_4param"),
    "7param":         (PARAM_NAMES_7,            "atom_counts_7param"),
    "hybrid":         (PARAM_NAMES_HYBRID,       "atom_counts_hybrid"),
    "bondorder":      (PARAM_NAMES_BONDORDER,    "atom_counts_bondorder"),
    "bondorder_ext":  (PARAM_NAMES_BONDORDER_EXT,"atom_counts_bondorder_ext"),
    "bondorder_ar":   (PARAM_NAMES_BONDORDER_AR, "atom_counts_bondorder_ar"),
    "extended":       (PARAM_NAMES_EXTENDED,     "atom_counts_extended"),
    "neighbour":      (PARAM_NAMES_NEIGHBOUR,    "atom_counts_neighbour"),
}
