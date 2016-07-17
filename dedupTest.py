# This program tests a cost function related to VM deduplication.

from math import ceil
from math import log
from random import randint

# number of physical machines
p = 1000
# number of failed machines
d = 10
# number of unique data segments
N = 100000
# max number of replicas among N segments
max_rep = 10000
# baseline replica No.
base_r = 5

## New Parameters:
# weight on survival
alpha = 0.5

# number of segments (need to be initialized)
n_seg = 0

'''
Initialize the list storing the numbers of replications for N segments.
We use Zipf distribution.
'''
def init_rep_list():
    global n_seg, max_rep, N
    l = []
    for i in range(N):
        n_i = ceil(max_rep/(i+1))
        l.append(n_i)
        n_seg += n_i
    return l

rep_list = init_rep_list()

'''
Probability of survival.
'''
def A(r):
    global d,p
    if r <=0:
        print('Invalid replica equals to 0.')
        exit(2)
    if r>d:
        return 1.0
    pr = 1.0
    for i in range(r):
        pr = pr*(d-i)/(p-i)
    return 1.0-pr

'''
Survival
'''
def svv(r):
    global N,alpha,base_r,n_seg
    if len(r) != N:
        print('Strategy profile length error.')
        exit(1)
    prod_Ar = 1.0
    for i in range(N):
        prod_Ar *= A(r[i])

    return prod_Ar

'''
Dedup rate
'''
def ddp(r):
    global N,alpha,base_r,n_seg
    if len(r) != N:
        print('Strategy profile length error.')
        exit(1)
    sum_r = 0
    for i in range(N):
        sum_r += r[i]

    return 1-sum_r/n_seg/base_r

'''
Cost function.
'''
def cost(r):
    global N,alpha,base_r,n_seg
    if len(r) != N:
        print('Strategy profile length error.')
        exit(1)
    sum_r = 0
    prod_Ar = 1.0
    for i in range(N):
        sum_r += r[i]
        prod_Ar *= A(r[i])
    if prod_Ar < 1e-7 or sum_r >= n_seg*base_r:
        return float('-inf')

    return alpha*log(prod_Ar) + (1-alpha)*log(1-sum_r/n_seg/base_r)

'''
Fluctuation.
'''
def fluc(r):
    global rep_list,base_r
    new_r = []
    for i in range(len(r)):
        nr = r[i]+randint(-1,1)
        if nr<base_r:
            new_r.append(base_r)
        elif nr>rep_list[i]*base_r:
            new_r.append(rep_list[i]*base_r)
        else:
            new_r.append(nr)
    return new_r

'''
Test
'''
def test():
    global base_r,rep_list
    r = []
    # sub-optimal
    for i in range(N):
        if i<0.02*N:
            r.append(base_r)
        else:
            r.append(base_r*rep_list[i])
    sub_opt = cost(r)
    opt_r = r
    print('Sub-optimal:')
    print('Survival = ', svv(opt_r), '\t Dedup rate = ', ddp(opt_r))
    unchanged_times = 0
    while unchanged_times < 10:
        r = fluc(opt_r)
        c = cost(r)
        if c > sub_opt:
            sub_opt = c
            opt_r = r
            unchanged_times = 0
        else:
            unchanged_times += 1
    print('Heuristic opt: ')
    print('Survival = ', svv(opt_r), '\t Dedup rate = ', ddp(opt_r))
    f = open('res.txt', 'w')
    for i in opt_r:
        f.write(str(i)+'\n')
    f.close()
