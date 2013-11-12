sum_formula = "YMn1.75Fe0.25O5"
density = 2.85
SpaceGroup = 55


def get_cs(x=0.25, delta=1):
    import pyasf
    sp = pyasf.sp
    cs = pyasf.unit_cell(SpaceGroup)
    cs.add_atom("Y", (sp.S(0.1344), sp.S(0.1687), 0), 1)
    
    cs.add_atom("FeP", (sp.S(0.3898), sp.S(0.3534), sp.S("1/2")), 0, occupancy =      x*delta) # 4h pyramidal position
    cs.add_atom("MnP", (sp.S(0.3898), sp.S(0.3534), sp.S("1/2")), 0, occupancy =  1 - x*delta) # 4h pyramidal position
    
    cs.add_atom("FeO", (0, sp.S("1/2"), sp.S(0.251)), 0, occupancy =     x * (1 - delta)) # 4f octahedral position
    cs.add_atom("MnO", (0, sp.S("1/2"), sp.S(0.251)), 0, occupancy = 1 - x * (1 - delta)) # 4f octahedral position
    
    cs.add_atom("O1", (0, 0, sp.S(0.260)), 1)
    cs.add_atom("O2", (sp.S(0.1601), sp.S(0.4417), 0), 1)
    cs.add_atom("O3", (sp.S(0.1523), sp.S(0.4291), sp.S("1/2")), 1)
    cs.add_atom("O4", (sp.S(0.3930), sp.S(0.2028), sp.S(0.2415)), 1)
    
    cs.get_tensor_symmetry()
    cs.build_unit_cell()
    
    cs.subs[cs.a]=7.3121
    cs.subs[cs.b]=8.5397
    cs.subs[cs.c]=5.7115
    
    return cs