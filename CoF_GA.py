#coding=utf-8 
'''This is a genetic algorithm for Compute-and-Forward scheme'''
'''We ever try to utilize the DEAP software package'''
'''仿真信道为2×2信道，能量上限暂定100'''
'''
1、初始化，产生初始个体和种群，确定交叉、变异和进化代数等数值
2、计算每一个体适应度
3、选择操作
4、交叉与变异操作
5、满足循环停止条件则退出（？）
'''
import sys
sys.path.append("/usr/local/lib/python2.7/dist-packages")
#print sys.path
import random
import math
import time
from copy import copy
from copy import deepcopy
from operator import itemgetter

from sage.all import *
from random import *
from numpy import *
from CoF_basic import *
from CoF_LLL import CoF_compute_fixed_pow_flex#评估函数

#指定迭代代数NGEN，变异概率MUPB，能量上限P_con
#每个维度染色体数目mount，
NGEN=20
MUPB=0.05
P_con=4095
mount=20
chromolength=int(math.ceil(math.log(P_con+1,2)))
#产生随机信道矩阵
if is_set_H == True:
    H_a = set_H_a
else:
    set_random_seed() # to avoid producing the same H_a in different threads
    H_a = matrix.random(RR, M, L, distribution=RealDistribution('gaussian', 1))
#H_a = matrix(RR, 2, 2, [[0.979236523248108, -0.129396925980777], [0.594475351529458, 0.666023537533719]])
#十进制转二进制数序列
#输入为十进制数，输出为转换后二进制0,1序列
'''
def int2bin(pop):
    popbin=[]#输出数组，存储二进制数的0,1序列
    str=bin(pop)[2:]
    for i in range(0,10-len(str)):
        popbin.append(0)
    for i in range(0,len(str)):
        popbin.append(int(str[i]))
    return popbin
#二进制序列转十进制
#输入为（十位长）二进制序列，输出为十进制数
def bin2int(popbin):
    pop=0
    for i in range(0,len(popbin)):
        pop+=popbin[i]*pow(2, 9-i)
    return pop
'''
def int2bin(pop):
    popbin=[]
    str=bin(pop)[2:]
    for i in range(0,chromolength-len(str)):
        popbin.append(0)
    for i in range(0,len(str)):
        popbin.append(int(str[i]))
    return popbin

def bin2int(popbin):
    pop=0
    for i in range(0,len(popbin)):
        pop+=popbin[i]*pow(2,chromolength-i-1)
    return pop

#根据能量上限P_con计算精度为1时的染色体长度
#返回全局变量：染色体长度
'''
def chromolength(P_con):
    return math.ceil(math.log(P_con+1,2))
'''
#产生均匀分布的初始种群
#P_con暂设定为100，染色体长度设定为10
'''
def initialpop():
    Population=[]#存储初始种群
    for i in arange(1,1023,100):
        for j in arange(1,1023,100):
            pop=[int2bin(i),int2bin(j)]
            #print pop
            #print bin2int(pop[0])
            Population.append(pop)    
    return Population
'''
#种群产生函数，种群个体采用均匀分布
#P_con为种群个体染色体数值上限（能量上限）
#chromolength为染色体长度（二进制位数）
#mount为每个维度染色体数目
def initialpop(P_con,L,mount):
    Population=[]#存储初始种群
    #chromolength=math.ceil(math.log(P_con+1,2))#向上取整
    for i in range(int(pow(mount,L))):
        pop=[]
        for j in range(L):
            l=int(uniform(0,P_con))
            pop.append(int2bin(l))
            #print pop
            #print bin2int(pop[0])
        Population.append(pop)    
    return Population

#print initialpop(100, 2, 10, 20)

#评估函数
#仅仅返回最大支持速率
#输入Power
def evaluate(P_t):
     return CoF_compute_fixed_pow_flex(P_t, P_con, False, H_a,False)
'''
print evaluate((1023,150))
H_a=set_H_a
print evaluate((bin2int(Population[0][0]),bin2int(Population[0][1])))
'''
#计算个体适应度函数
def fitvalue(pop):
    fitness=[]
    for i in range(0,mount*mount):
        popvalue=[]
        for j in range(L):
            popvalue.append(bin2int(pop[i][j]))
        fitness.append(evaluate(popvalue))
    return fitness

#定义交叉函数
def crossover(ind1,ind2):
    #for i in range(2):
    #chromolength=math.ceil(math.log(P_con+1,2))
    p1=randint(0,chromolength)
    p2=randint(0,chromolength-1)
    if p2>=p1:
        p2+=1
    else:
        p1,p2=p2,p1
    p3=randint(0,chromolength)
    p4=randint(0,chromolength-1)
    if p4>=p3:
        p4+=1
    else:
        p3,p4=p4,p3        
    #ind1[i][p1:p2],ind2[i][p1:p2]=ind2[i][p1:p2],ind1[i][p1:p2]
    if L==2:
        return (ind1[0][:p1]+ind2[0][p1:p2]+ind1[0][p2:],ind1[1][:p3]+ind2[1][p3:p4]+ind1[1][p4:])\
            ,(ind2[0][:p1]+ind1[0][p1:p2]+ind2[0][p2:],ind2[1][:p3]+ind1[1][p3:p4]+ind2[1][p4:])
    elif L==3:
        p5=randint(0,chromolength)
        p6=randint(0,chromolength-1)
        if p6>=p5:
            p6+=1
        else:
            p5,p6=p6,p5
        return (ind1[0][:p1]+ind2[0][p1:p2]+ind1[0][p2:],ind1[1][:p3]+ind2[1][p3:p4]+ind1[1][p4:],ind1[2][:p5]+ind2[2][p5:p6]+ind1[2][p6:])\
            ,(ind2[0][:p1]+ind1[0][p1:p2]+ind2[0][p2:],ind2[1][:p3]+ind1[1][p3:p4]+ind2[1][p4:],ind2[2][:p5]+ind1[2][p5:p6]+ind2[2][p6:])
    else:
        raise Exception("error: cannot support  a channel more than 3 input &output")
    #return ind1,ind2
    
'''
#定义交叉函数
def crossover(ind1,ind2,indpb):
    ind=ind1
    for i in range(10):
        if random()<indpb:
            #ind1[0][i],ind2[0][i]=ind2[0][i],ind1[0][i]
            ind1[0][i]=ind2[0][i]
            ind2[0][i]=ind[0][i]
        if random()<indpb:            
            #ind1[1][i],ind2[1][i]=ind2[1][i],ind1[1][i]
            ind1[1][i]=ind2[1][i]
            ind2[1][i]=ind[1][i]
    #return ind1,ind2
'''
#定义变异函数，ind为个体，indpb为变异概率
#@parallel(ncpus=Cores)
def mutate(individual,indpb):
    for i in range(0,L):
        for j in range(chromolength):
            if random() < indpb:
                individual[i][j] = type(individual[i][j])(not individual[i][j])
    #return individual

#排序函数，按照适应度对种群进行排序
#返回排序后的种群序列(zip,sorted)
#函数中对fitscore乘以100,是为了提高各个体之间适应度差异
def rankPop(pop):
    fitness=fitvalue(pop)
    sum_fit=sum(fitness)
    fitscore=[]
    for i in range(len(fitness)):
        fitscore.append(100*fitness[i]/sum_fit)
    pairedpop=zip(pop,fitscore)
    rankedpop=sorted(pairedpop,key=itemgetter(-1),reverse=True)
    return rankedpop


#使用轮盘赌选择法
#输入归一化的适应度
#输出相应选择项的下标
def Roulette(fitscore):
    caculatefitscore=0
    for i in range(len(fitscore)):
        caculatefitscore+=fitscore[i]
        if caculatefitscore>100*random():
            return i
        
#产生新种群函数
#采用精英选择和轮盘赌的方法
#输入当代种群
#输出新一代种群
def NewPop(pop):
    lenpop=len(pop)
    rankedpop=rankPop(pop)
    fitscore=[item[-1] for item in rankedpop]
    rankedchromes=[item[0] for item in rankedpop]    
    newpop=[]
    #精英选择，选择当代适应度最高的前百分之十
    newpop.extend(rankedchromes[0:lenpop/10])
    #轮盘赌选择
    #选择的个体参与遗传操作：交叉与变异
    while len(newpop)!=lenpop:
        ind1,ind2=[],[]
        index1=Roulette(fitscore)
        index2=Roulette(fitscore)
        while index2==index1:
            index2=Roulette(fitscore)
        ind1=rankedchromes[index1]
        ind2=rankedchromes[index2]
        ind1,ind2=crossover(ind1, ind2)
        #crossover(ind1, ind2, CXPB)
        #ind1=mutate(ind1,MUPB)
        #ind2=mutate(ind2,MUPB)
        mutate(ind1, MUPB)
        mutate(ind2, MUPB)
        newpop.extend([ind1,ind2])
    return newpop
    
#@parallel(ncpus=Cores)
def GeneticAlgorithm():
    #产生初始种群
    pop=initialpop(P_con,L,mount)
    for i in range(0,NGEN):
        fitness=fitvalue(pop)
        print "Generation: ", i
        print "该代最大适应度：",max(fitness)
        print "最大适应度对应个体：",pop[fitness.index(max(fitness))]

        #选择并交叉变异产生下一代种群
        newpop=[]
        newpop=NewPop(pop)
        pop=copy(newpop)
    fitness=fitvalue(pop)
    print "最大适应度：",max(fitness)
    print"最大适应度对应个体：",pop[fitness.index(max(fitness))],

if __name__=="__main__":
    t1=time.time()
    GeneticAlgorithm()
    t2=time.time()
    print "程序耗时：%i 秒"%(t2-t1)