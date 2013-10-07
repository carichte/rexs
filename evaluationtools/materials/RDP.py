sum_formula = "H2O4P1Rb1"
density = 2.85
SpaceGroup = 122

def get_cs():
    import StructureTensor
    sp = StructureTensor.sp
    cs = StructureTensor.unit_cell(SpaceGroup)
    cs.add_atom("Rb", (0,0,sp.S("1/2")), 0)
    cs.add_atom("P", (0, 0, 0), 1)
    cs.add_atom("O", (sp.S(0.1421), sp.S(0.0863), sp.S(0.1199)), 1)
    cs.add_atom("H", (sp.S(0.1382), sp.S(0.223) ,sp.S("1/8")), 1)
    
    
    cs.get_tensor_symmetry()
    cs.build_unit_cell()
    #cs.subs[cs.a]=7.622
    cs.subs[cs.a]=7.6247 # Ovch
    #cs.subs[cs.a]=7.6065 # Med T
    #cs.subs[cs.a]=7.591 # Low T
    #cs.subs[cs.a]=7.591 # Low T
    #cs.subs[cs.c]=7.315
    cs.subs[cs.c]=7.3041 # Ovch
    #cs.subs[cs.c]=7.286 # Med T
    #cs.subs[cs.c]=7.257 # Low T
    #cs.subs[cs.c]=7.257 # Low T
    #cs.subs["x"]=0.30499
    return cs