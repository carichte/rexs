sum_formula = "TiO2"
SpaceGroup = 136

def get_cs():
    import pyasf
    sp = pyasf.sp
    x = sp.Symbol("x", real=True, unbounded=False)
    cs = pyasf.unit_cell(SpaceGroup)
    cs.add_atom("Ti", (0,0,0), 0)
    cs.add_atom("O", (x, x, 0), 1)
    
    
    cs.subs[cs.a]=4.5924
    cs.subs[cs.c]=2.9575
    cs.subs[x]=0.30499
    cs.get_tensor_symmetry()
    cs.build_unit_cell()
    return cs