import jpype
import jpype.imports
import os

# Set up variables
jar_path = "gtp.jar"  
newick1 = "(A:1,B:1,C:1);"
newick2 = "(B:1,A:1,C:1);"
rooted = True

# Start JVM
jpype.startJVM(classpath=[jar_path])

# Import classes
PolyMain = jpype.JClass("polyAlg.PolyMain")
PhyloTree = jpype.JClass("distanceAlg1.PhyloTree")

# Create PhyloTree objects
tree1 = PhyloTree(newick1, rooted)
tree2 = PhyloTree(newick2, rooted)

# Compute distance
outFile = None  # No output file needed
distance = PolyMain.getGeodesic(tree1, tree2, outFile).getDist()

print("BHV distance:", distance)

# Shutdown JVM
jpype.shutdownJVM()
