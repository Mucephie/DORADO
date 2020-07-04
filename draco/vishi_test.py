num = int(input("give a num: "))
numrange = range(2,num)

for i in numrange:
        k=0
        for j in range(2,num-1):
                if ((i%j)==0):
                        # k=0
                        break
                else:
                        k += 1
        if (k>0):
                print(i, " is a prime number")
else:
        print("over")