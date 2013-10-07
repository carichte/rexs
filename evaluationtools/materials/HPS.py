import os
import numpy as np 

sum_formula = "Ho2PdSi3"
ciffile = "HPS_D1.cif"
ciffile = os.path.join(os.path.split(__file__)[0], ciffile)
#density = 4.63
SpaceGroup = 1
layertypes = ("a", "b", "c", "d")


def add_layer(cs, label, zoffset=0):
    if label.lower() not in layertypes:
        raise ValueError("Layer with label %s not supported"%label)
    fname = os.path.join(os.path.split(__file__)[0], "%s.txt"%label.upper())
    layer = np.loadtxt(fname, dtype=str)
    for atom in layer:
        label, x, y, z = atom
        z = str(float(z)/8 + zoffset)
        label += str(cs.elements.values().count(label))
        cs.add_atom(label, (x,y,z))


def get_cs(stack = None):
    import StructureTensor
    sp = StructureTensor.sp
    if stack==None:
        cs = StructureTensor.unit_cell(ciffile)
    else:
        cs = StructureTensor.unit_cell(1)
        for i in range(len(stack)):
            add_layer(cs, stack[i], float(i)/len(stack))
        
        cs.subs[cs.a] = 8.1
        cs.subs[cs.b] = 8.1
        cs.subs[cs.c] = 32.0/8*len(stack)
        cs.subs[cs.alpha] = sp.pi/2
        cs.subs[cs.beta] = sp.pi/2
        cs.subs[cs.gamma] = 2*sp.pi/3
    
    cs.build_unit_cell()
    return cs