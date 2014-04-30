import os
sum_formula = "H2O4P1Rb1"
density = 2.85
SpaceGroup = 122
SpaceGroupLT = 43
ciffile = "cif/Fe3Ga_Fm3m.cif"

ciffile = os.path.join(os.path.split(__file__)[0], ciffile)

def get_cs():
    import pyasf
    sp = pyasf.sp
    cs = pyasf.unit_cell(ciffile, resonant="Fe")
    cs.get_tensor_symmetry()
    cs.build_unit_cell()
    return cs