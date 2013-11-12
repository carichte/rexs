sum_formula = "BaTiO3"
density = 5.85
SpaceGroup = 99


def get_cs():
    import pyasf
    sp = pyasf.sp
    np = pyasf.np
    
    delta_1 = sp.Symbol("delta_Ti", real=True)
    delta_2 = sp.Symbol("delta_O1", real=True)
    delta_3 = sp.Symbol("delta_O2", real=True)

    cs = pyasf.unit_cell(SpaceGroup)
    
    cs.add_atom("Ba", (0,0,0), 1, dE=-4)
    cs.add_atom("Ti", (sp.S("1/2"),sp.S("1/2"),sp.S("1/2") - delta_1), 1, dE=14.5)
    cs.add_atom("O1", (sp.S("1/2"),sp.S("1/2"),0 + delta_2), 1)
    cs.add_atom("O2", (sp.S("1/2"),0,sp.S("1/2") + delta_3), 1)
    
    cs.subs[cs.a] = 3.9998
    cs.subs[cs.c] = 4.0180
    
    cs.subs[delta_1] = 0.018
    cs.subs[delta_2] = 0.016
    cs.subs[delta_3] = 0.015
    
    cs.get_tensor_symmetry()
    cs.build_unit_cell()
    return cs