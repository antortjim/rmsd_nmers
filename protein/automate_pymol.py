from pymol import cmd

# Automation of simple Pymol tasks
# Takes a picture of the alignment between 2 fragments
# One fragment coming from s1 (structure 1) in c1 (chain in s1), between residues i1 and i2
# The other coming from s2 (structure 2) in c2 (chain in s2), between residues j1 and j2

## BUG: Does not remove fragments between the same resi but in chains different from c1/c2
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
