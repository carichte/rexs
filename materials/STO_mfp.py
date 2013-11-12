sum_formula = "SrTiO3"
density = 5.13
SpaceGroup = 99


def get_cs():
    import pyasf
    sp = pyasf.sp
    np = pyasf.np
    
    delta_1 = 0
    delta_2 = 0
    delta_3 = 0

    cs = pyasf.unit_cell(SpaceGroup)
    
    cs.add_atom("Sr", (0,0,0), 1, dE=-4)
    cs.add_atom("Ti", (sp.S("1/2"),sp.S("1/2"),sp.S("1/2") - delta_1), 1, dE=14.5)
    cs.add_atom("O1", (sp.S("1/2"),sp.S("1/2"),0 + delta_2), 1)
    cs.add_atom("O2", (sp.S("1/2"),0,sp.S("1/2") + delta_3), 1)
    
    cs.subs[cs.a] = 3.905 - 1.3e-4
    cs.subs[cs.c] = 3.905 + 5e-3
    
    cs.subs[delta_1] = 0.0
    cs.subs[delta_2] = 0.0
    cs.subs[delta_3] = 0.0
    
    cs.get_tensor_symmetry()
    cs.build_unit_cell()
    return cs