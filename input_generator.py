import sys
import random

SIZE = int(sys.argv[1])

print(str(SIZE) + " " + str(SIZE - 1))
for i in range(0, SIZE):
    print (str(random.randint(100,1000)) + " " + str(random.randint(10,100)))