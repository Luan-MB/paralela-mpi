import sys
import random

profitList = []
weigthList = []

SIZE = int(sys.argv[1])
for i in range(0, SIZE):
    n = random.randint(10,100)
    m = random.randint(10,100)
    profitList.append(n)
    weigthList.append(m)

print(str(SIZE) + " " + str(random.randint(1000, 5000)))
for i in range(0, SIZE):
    print (str(profitList[i]) + " " + str(weigthList[i]))