import time
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool 

from individual import individual 
from config_var import config



class generation:

    get = config().get

    __neighbourhood = get("neighbourhood") # __neighbourhood = get("neighbourhood") #规则同上一层的左右三个细胞有关
    __gene_length = 2**(2*__neighbourhood+1) 

    __neighbour_radius = get("neighbour_radius")  #邻居域半径
    __neighbour_number = (2*__neighbour_radius+1)*(2*__neighbour_radius+1)
    __group_width = get("group_width")
    __group_height = get("group_height")
    __group_number = __group_width*__group_height
    # __parent_number = __group_number
    # __child_number = __group_number

    __generation_num = get("generation_num")  #运行多少代

    __pool_size = (__group_number,__gene_length)#基因池大小
    __parent_gene_pool = np.zeros(__pool_size,dtype=np.int) #当前种群基因池
    __child_gene_pool = np.zeros(__pool_size,dtype=np.int)
 
    # __selected_parents = np.zeros(__child_number*2) #选出的父母

    __parent_adapt = np.zeros(__group_number)

    __adapt_max_index = 0  #适应度最佳的个体的标号
    __adapt_max = []


    def __init__(self):
        #初代，随机初始化基因池
        #写到文件里
        
        initial_flag = np.loadtxt('config')

        if initial_flag:      #基因池进行随机初始化
            print("————————————————————————————————")
            print("initial_flag为1,基因池进行随机初始化")
            gene_pool_initial = np.random.randint(2,size =(self.__pool_size))
            np.savetxt('gene_pool',gene_pool_initial,fmt='%d')
            np.savetxt('config',[0],fmt='%d')
            np.savetxt('adapt',[],fmt='%d')

        else:                 #基因池不再初始化，直接读取
            print("————————————————————————————————")
            print("initial_flag为0,基因池不再初始化")
            #读取上回训练的基因池
            gene_pool_initial = np.loadtxt('gene_pool')
        
        self.__parent_gene_pool = gene_pool_initial.astype(int)
        # self.__parent_gene_pool = np.random.randint(2,size =(self.__pool_size))
        
        # for i in range(self.__generation_num):
        #     self.select_parents()
        #     self.give_birth()
        #     self.__parent_gene_pool =  self.__child_gene_pool#当前种群基因池
        #     # self.__child_gene_pool = np.random.randint(2,size =(self.__pool_size))
        
        # #基因池写到文件里，结束该世代
        # np.savetxt('gene_pool',self.__parent_gene_pool)
        
    def get_neighbours(self,k):
        width = self.__group_width
        height = self.__group_height
        radius = self.__neighbour_radius
        #一维转二维
        i = int(k/width)
        j = k % width

        neighbours = []
        #二维生境的首尾相连
        for x in range(-radius,radius+1):
            if i+x < 0 or i+x >= height:
                i_x = (i+x+height)%height  
            else:
                i_x = i+x
            for y in range(-radius,radius+1):
                if j+y < 0 or j+y >= width:
                    j_y = (j+y+width)%width
                else:
                    j_y = j+y
                neighbours.append((i_x)*width+j_y)
        # print(k)
        # print(neighbours)
        return neighbours

    def get_adaptability(self,k):
        indiv = individual(self.__parent_gene_pool[k])
        self.__parent_adapt[k] = indiv.adaptability
        print("get_adaptability",k,":",self.__parent_adapt[k])
        print("get_adaptability",k,":",self.__parent_adapt)


    def give_birth(self,k):
        print("give_birth",k)

        parents = self.select_parents(k)
        child_gene = self.crossover_and_mutation(parents)
        
        child = individual(child_gene)
        child_adapt = child.adaptability
        #如果孩子的适应度大于父母才替换，否则保持
        if child_adapt > self.__parent_adapt[k]:
            self.__child_gene_pool[k]=child_gene
        else:
            self.__child_gene_pool[k]=self.__parent_gene_pool[k]

        

    #通过父代适应度产生父母轮盘，用轮盘赌产生选择的父母
    #轮盘赌的实现方式：记录一个数列，
    # [adapt1,adapt1+2,adapt1+2+3...adapt_sum] 
    # random[parent_number]，落到第i个区间就选择第i个做父母，共选择parent_number 个父母，
    def select_parents(self,k):

        adapt_roulette = []
        neighbours = self.get_neighbours(k)
        neighbour_number = self.__neighbour_number
        parents = []
        
        adapt_sum = 0
        for n in neighbours:
            adapt = self.__parent_adapt[n]

            adapt_sum = adapt_sum + adapt
            adapt_roulette.append(adapt_sum)

        #从邻居中选两个父母
        rand = np.random.uniform(0,adapt_sum,2)
        for i in range(2):
            for j in range(neighbour_number):
                    if rand[i] < adapt_roulette[j]:
                        parents.append(neighbours[j])
                        # print("i=",i," selected parents:",self.__selected_parents[i])
                        break
        return parents

    #一点交叉，把父母基因赋予一个子代，并随机突变 
    def crossover_and_mutation(self,parents):
        point = np.random.randint(0,self.__gene_length)
        mather_gene = self.__parent_gene_pool[int(parents[0])]
        father_gene = self.__parent_gene_pool[int(parents[1])]

        child_gene = np.zeros(self.__gene_length,dtype=np.int)
        child_gene[0:point] = mather_gene[0:point]
        child_gene[point:self.__gene_length] = father_gene[point:self.__gene_length]
        mutation_num = np.random.randint(6)
        for j in range(mutation_num):
            mutation_point = np.random.randint(self.__gene_length)
            child_gene[mutation_point] = (child_gene[mutation_point]+1)%2  #  0，1互换
            
        return child_gene

    def save_gene(self):
        self.__parent_gene_pool =  self.__child_gene_pool#当前种群基因池
        np.savetxt('gene_pool',self.__parent_gene_pool)
        
    def save_adapt(self):
        print(self.__parent_adapt)
        adapt_max = 0
        for adapt in self.__parent_adapt:
            if adapt > adapt_max:
                adapt_max = adapt
        
        with open('adapt', 'a') as f:
        # f.write(adapt)
            np.savetxt(f,[adapt_max])
    
    def print_adapt(self):
        print(self.__parent_adapt)
    

if __name__ == '__main__':

    get = config().get
    group_width = get("group_width")
    group_height = get("group_height")
    group_number = group_width*group_height

    Gen1 = generation()
    


    stime = time.time()
    #python multiprocessing 并行计算
    #windows 下只能在main函数里写
    for step in range(1):
        k = range(group_number)
        pool = Pool()
        pool.map(Gen1.get_adaptability, k)
        pool.close()
        pool.join()

        Gen1.print_adapt()

        pool = Pool()
        pool.map(Gen1.give_birth, k)
        pool.close()
        pool.join()
        
        Gen1.save_adapt()
        Gen1.save_gene()
        
    # Gen1.save_data()
    
    etime = time.time() - stime
        
    print('并行用时:{}秒'.format(etime))

    stime = time.time()

    for step in range(1):
        for k in range(group_number):
            Gen1.get_adaptability(k)
        for k in range(group_number):    
            Gen1.give_birth(k)
        Gen1.save_adapt()
        Gen1.save_gene()

    etime = time.time() - stime
    print('串行用时:{}秒'.format(etime))
