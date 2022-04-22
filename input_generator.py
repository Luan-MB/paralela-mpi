import sys
import random

profitList = []
weigthList = []

SIZE = int(sys.argv[1])

print(str(SIZE) + " " + str(SIZE))
for i in range(0, SIZE):
    print (str(random.randint(100,1000)) + " " + str(random.randint(10,100)))