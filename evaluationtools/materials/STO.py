sum_formula = "SrTiO3"
density = 5.13
SpaceGroup = 221


def get_cs():
    import StructureTensor
    sp = StructureTensor.sp
    np = StructureTensor.np
    
    delta_1 = 0
    delta_2 = 0
    delta_3 = 0

    cs = StructureTensor.unit_cell(SpaceGroup)
    
    cs.add_atom("Sr", (0,0,0), 1, dE=-4)
    cs.add_atom("Ti", (sp.S("1/2"),sp.S("1/2"),sp.S("1/2") - delta_1), 1, dE=14.5)
    cs.add_atom("O", (sp.S("1/2"),sp.S("1/2"),0 + delta_2), 1)
    
    cs.subs[cs.a] = 3.905
    #cs.subs[cs.c] = 
    
    #cs.subs[delta_1] = 0.018
    #cs.subs[delta_2] = 0.016
    #cs.subs[delta_3] = 0.015
    
    cs.get_tensor_symmetry()
    cs.build_unit_cell()
    return cs