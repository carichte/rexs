import evaluationtools as et
sum_formula = "LiNbO3"
density = 4.63
SpaceGroup = 161

# Expansion Coefficients (Citation: J. Appl. Phys. 40, 4637 (1969); doi: 10.1063/1.1657244)
TExp = et.TempExpansion()
TExp.alpha = ((1.44e-5, 1.54e-5, 0), (1.54e-5, 1.59e-5, 0), (0, 0, 0.75e-5))
TExp.beta  = (( 7.1e-9,  5.3e-9, 0), (5.3e-9,   4.9e-9, 0), (0, 0, -7.7e-9))

def get_cs(**kwargs):
    import pyasf
    sp = pyasf.sp
    cs = pyasf.unit_cell(SpaceGroup)
    
    cs.add_atom("Li", (0, 0, 0.2829), 1)
    cs.add_atom("Nb", (0, 0, 0), 0)
    cs.add_atom("O", (0.0492, 0.3446, 0.0647), 1)
    cs.subs[cs.a] = 5.14829
    cs.subs[cs.c] = 13.8631
    
    if kwargs.has_key("temperature"):
        fac = TExp(kwargs["temperature"])
        cs.subs[cs.a] *= fac[0,0]
        cs.subs[cs.c] *= fac[2,2]
    
    
    cs.get_tensor_symmetry()
    cs.build_unit_cell()
    return cs