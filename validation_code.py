sequential=[]
parallel=[]
with open("sequential.txt") as fp: 
    while True: 
        #count += 1
        line = fp.readline() 
  
        if not line: 
            break
        #print(line)
        sequential.append(list(map(int,line.split())))
with open("parallel.txt") as fp: 
    while True: 
        #count += 1
        line = fp.readline() 
  
        if not line: 
            break
        #print(line)
        parallel.append(list(map(int,line.split())))
#print(sequential)
#print(parallel)
if(parallel==sequential):
	print("valid")
else:
	print("invalid")
#print(parallel==sequential)