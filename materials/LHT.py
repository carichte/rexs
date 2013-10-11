sum_formula = "NH2.5Ti1.87O4"
density = 1.367
SpaceGroup = 71

def get_cs():
    import StructureTensor
    sp = StructureTensor.sp
    cs = StructureTensor.unit_cell(SpaceGroup)
    cs.add_atom("Ti", (0,0.3048,sp.S("1/2")), 0)
    cs.add_atom("O1", (0,0.2208,0), 1)
    cs.add_atom("O2", (0,0.3619,0), 1)
    cs.add_atom("N",  (0,0, 0.230), 1)
    
    cs.get_tensor_symmetry()
    cs.build_unit_cell()
    
    cs.subs[cs.a]=3.7926
    cs.subs[cs.b]=18.458
    cs.subs[cs.c]=2.9774
    return cs