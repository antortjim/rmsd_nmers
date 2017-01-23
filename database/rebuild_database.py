import Bio.PDB as PDB
import time

pdbList = PDB.PDBList()
pdb_textfn = "../pdb.txt"


with open(pdb_textfn) as pdb_textfile:
    # PDB_textfn contains several lines where everyline is a PDB id + chain id
    # For example
    #   2ptcE
    #   2ptcI
    #   XXXXY
    for line in pdb_textfile:
       # Select PDB ID
       pdb_id = line[:4].lower()
       # Select chain
       #chain = line[4]
       print "Query %s" % pdb_id
       try:
           pdb_path = pdbList.retrieve_pdb_file(pdb_id, pdir = ".")
       
       except IOError:
               try:
                  time.sleep(3)
                  pdb_path = pdbList.retrieve_pdb_file(pdb_id,
                                                       obsolete = True,
                                                       pdir = ".")
               except IOError:
                  raise
       # Pass pdb_path (pdb file path as returned by retrieving function)
       # and chain (parsed from pdb.txt input file) to make_pdb()
       #splitter.make_pdb(pdb_path, chain)
     
                
    
