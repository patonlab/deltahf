"""Physical constants and model definitions shared across deltahf modules."""

HARTREE_TO_KCAL = 627.5094740631
HARTREE_TO_EV = 27.211386245988

# Atom type names for each classification scheme
PARAM_NAMES_4 = ["C", "H", "N", "O"]
PARAM_NAMES_7 = ["C", "H", "N", "O", "C_prime", "N_prime", "O_prime"]
PARAM_NAMES_HYBRID = ["C_sp3", "C_sp2", "C_sp", "H", "N_sp3", "N_sp2", "N_sp", "O_sp3", "O_sp2", "O_sp"]
PARAM_NAMES_EXTENDED = [
    "C_sp3_3H", "C_sp3_2H", "C_sp3_1H", "C_sp3_0H",
    "C_sp2_2H", "C_sp2_1H", "C_sp2_0H", "C_sp",
    "H", "N_sp3", "N_sp2", "N_sp", "O_sp3", "O_sp2", "O_sp",
]

# Maps model name -> (param_names, MoleculeResult atom_counts field)
MODEL_DEFS = {
    "4param": (PARAM_NAMES_4, "atom_counts_4param"),
    "7param": (PARAM_NAMES_7, "atom_counts_7param"),
    "hybrid": (PARAM_NAMES_HYBRID, "atom_counts_hybrid"),
    "extended": (PARAM_NAMES_EXTENDED, "atom_counts_extended"),
}
