import evaluationtools as et
sum_formula = "LiTaO3"
density = 7.41
# Expansion Coefficients (Citation: J. Appl. Phys. 40, 4637 (1969); doi: 10.1063/1.1657244)
TExp = et.TempExpansion()
TExp.alpha = ((1.61e-5, 1.62e-5, 0), (1.62e-5, 1.54e-5, 0), (0, 0, 0.22e-5))
TExp.beta  = (( 7.5e-9,  5.9e-9, 0), (5.9e-9,   7.0e-9, 0), (0, 0, -5.9e-9))

SpaceGroup = 161


def get_cs(**kwargs):
    import StructureTensor
    sp = StructureTensor.sp
    np = StructureTensor.np
    cs = StructureTensor.unit_cell(SpaceGroup)
    
    cs.add_atom("Li", (0, 0, 0.2821), 1, assume_complex=True)
    cs.add_atom("Ta", (0, 0, 0), 0, assume_complex=True)
    cs.add_atom("O",  (0.0534, 0.3396, 0.0695), 1, assume_complex=True)
    cs.subs[cs.a] =  5.15428
    cs.subs[cs.c] = 13.7835
    if kwargs.has_key("temperature"):
        fac = TExp(kwargs["temperature"])
        cs.subs[cs.a] *= fac[0,0]
        cs.subs[cs.c] *= fac[2,2]
    
    cs.get_tensor_symmetry()
    cs.build_unit_cell()
    
    return cs