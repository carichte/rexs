sum_formula = "H2O4P1Rb1"
density = 2.85
SpaceGroup = 136

def get_cs():
    import StructureTensor
    sp = StructureTensor.sp
    x = sp.Symbol("x", real=True, unbounded=False)
    cs = StructureTensor.unit_cell(SpaceGroup)
    cs.add_atom("Ti", (0,0,0), 0)
    cs.add_atom("O", (x, x, 0), 1)
    
    
    #cs.get_tensor_symmetry()
    cs.build_unit_cell()
    cs.subs[cs.a]=4.5924
    cs.subs[cs.c]=2.9575
    cs.subs[x]=0.30499
    return cs