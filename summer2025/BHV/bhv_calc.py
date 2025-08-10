import jpype

# Start JVM with the JAR
jpype.startJVM(classpath=["gtp.jar"])

# Import Java classes
PhyloTree = jpype.JClass("distanceAlg1.PhyloTree")
PolyMain = jpype.JClass("polyAlg.PolyMain")

def calculate_bhv_distance(newick1, newick2, rooted=True):
    """Calculate BHV distance between two Newick format trees."""
    tree1 = PhyloTree(newick1, rooted)
    tree2 = PhyloTree(newick2, rooted)
    return PolyMain.getGeodesic(tree1, tree2, None).getDist()

# Example usage
if __name__ == "__main__":
    distance = calculate_bhv_distance("((A,B),C);", "((A,C),B);")
    print(f"BHV Distance: {distance}")
    jpype.shutdownJVM()