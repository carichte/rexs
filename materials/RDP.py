import os
sum_formula = "H2O4P1Rb1"
density = 2.85
SpaceGroup = 122
SpaceGroupLT = 43

def get_cs(T = 298., LT = False):
    import pyasf
    sp = pyasf.sp
    if LT:
        ciffile = "cif/RDP_orthorhombic_146K.cif"
        ciffile = os.path.join(os.path.split(__file__)[0], ciffile)
        cs = pyasf.unit_cell(ciffile, resonant="Rb")
    else:
        ciffile = "cif/RDP_tetragonal_RT.cif"
        ciffile = os.path.join(os.path.split(__file__)[0], ciffile)
        cs = pyasf.unit_cell(ciffile, resonant="Rb")
        alpha_x = 0.0002262773722627715 # expansion coefficient
        alpha_z = 0.0004233576642335819 # expansion coefficient
        cs.subs[cs.a] +=  alpha_x * (T - 298.)
        cs.subs[cs.c] +=  alpha_z * (T - 298.)
    cs.get_tensor_symmetry()
    cs.build_unit_cell()
    return cs