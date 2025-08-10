import jpype

# Start JVM
jpype.startJVM(classpath=["test.jar"])  # use gtp.jar later

Test = jpype.JClass("Test")
result = Test.add(2, 3)
print("Test add result:", result)

jpype.shutdownJVM()
