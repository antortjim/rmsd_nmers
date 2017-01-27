from pymol import cmd

def show_fragment(s1, s2, c1, i1, i2, c2, j1, j2):
    cmd.load(s1 + ".ent")
    cmd.load(s2 + ".ent")


    cmd.remove("solvent")

    cmd.remove("model %s and not resi %d:%d" % (s1, i1, i2))
    cmd.remove("model %s and not resi %d:%d" % (s2, j1, j2))

    cmd.remove("model %s and not chain %s" % (s2, c1))
    cmd.remove("model %s and not chain %s" % (s2, c2))

    cmd.align(s1, s2, cutoff=2.0,
            cycles=5, gap=-10.0,  extend=-0.5,
            max_gap=50, object=None, matrix='BLOSUM62',
            mobile_state=0,  target_state=0,  quiet=1,
            max_skip=0,  transform=1,  reset=0)
    
    
    cmd.bg_color("white")

    cmd.hide("lines")
    cmd.show("cartoon")

  
    cmd.zoom(s1)
    
    cmd.png("../plots/sample.png")


cmd.extend('show_fragment', show_fragment)
