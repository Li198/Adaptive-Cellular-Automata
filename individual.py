import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl 
from multiprocessing import Pool 

from config_var import config


class individual:
    #定义基本属性

    get = config().get

    __neighbourhood = get("neighbourhood") # __neighbourhood = get("neighbourhood") #规则同上一层的左右三个细胞有关
    __width = get("width")
    __height = get("height")
    __test_step =  get("individual_test_step")   #通过多少步来得到 adaptability

    __gene = []

    __cells = np.zeros((__height,__width),dtype=np.int)

    adaptability = 0



    #定义构造方法

    def __init__(self,gene):
        self.__gene = gene
        #计算 测试 的时间

        self.__test_adaptability()
      
        
        
        
        
        


    #跑 __height 时间步
    def __run_cells(self):
        for i in range(1,self.__height):
            for j in range(0,self.__width):
                self.__cells[i][j] = self.__use_gene(i,j)


        # width = self.__width
        # height = self.__height

        #计算第 1 到最后一层的cell
        # for i in range(1,height):
        #     k = [i*width+j for j in range(width)]
        #     # print(k)
        #     use_gene = self.__use_gene
        #     pool = Pool()
        #     pool.map(use_gene, k)
        #     pool.close()
        #     pool.join()




    #判别方法，目前是上一层的 共7各相邻细胞
    def __use_gene(self,i,j):
        # i = int(k/self.__width)
        # j = k % self.__width

        pre_layer = self.__cells[i-1]
        gene_index = 0
        #  (j+k+self.__width)%self.__width   循环，首尾相连
        #  gene_index是把cell的邻域转成二进制串，再转成十进制后的基因坐标
        for k in range(-self.__neighbourhood,self.__neighbourhood-1,1 ):
            gene_index = gene_index*2 + pre_layer[(j+k+self.__width)%self.__width] 
        return self.__gene[int(gene_index)]


    def __test_gene(self):
        start = self.__cells[0]
        end = self.__cells[-1]
        threshold = (self.__width-1)/2
        #开始和结束时，0多还是1多
        winner_start = 1 if np.count_nonzero(start) > threshold else 0    
        winner_end = 1 if np.count_nonzero(end) > threshold else 0
        #结束时，1多就数1
        winner_propertion = np.count_nonzero(end) / self.__width
        if winner_end == 0:
            winner_propertion = 1-winner_propertion
        #初始占多数的细胞最终占多数,返回所占比例
        if winner_start == winner_end :   
            # print("winner_propertion",winner_propertion)
            return winner_propertion
        else:     #初始占多数的细胞最终占少数
            # print("winner_propertion",0)
            return -winner_propertion

    def __test_adaptability(self):
        sum = 0
        step = self.__test_step

        for i in range(step):
            self.__cells[0] = np.random.randint(2,size=self.__width)
            self.__run_cells()

            # print("test",i)
            sum  = sum + self.__test_gene()
        self.adaptability = (sum/self.__test_step) + 1     #适应性为 0—2 的一个值，越大越好

    def show_cells(self):

        self.__cells[0] = np.random.randint(2,size=self.__width)



        __cells = self.__cells
        __adapt = self.adaptability
        title = "adapt:"+str(__adapt)
        plt.title(title)
        plt.rcParams['figure.figsize'] = (6.0, 6.0)# 设置figure_size尺寸
        # colors = ['white', 'black'] 
        # cmap = mpl.colors.ListedColormap(colors)
        # plt.imshow(__cells,cmap=cmap)
        plt.imshow(__cells)

        plt.show()

    def test(self):

        self.__cells[0] = np.random.randint(2,size=self.__width)
        
        self.__run_cells()

    
    
# def final():
#     gene = [0,0,0,0,0,0,0,0,0,0,0,0,
#     0,0,0,0,0,0,0,0,0,0,0,0,
#     0,0,0,0,0,0,0,0,0,0,0,0,
#     0,0,0,0,0,0,0,0,0,0,0,0,
#     0,0,0,0,0,0,0,0,0,0,0,0,
#     0,0,0,0,0,0,0,0,0,0,0,0,
#     0,0,0,0,0,0,0,0,0,0,0,0,
#     0,0,0,0,0,0,0,0,0,0,0,0,
#     0,0,0,0,0,0,0,0,0,0,0,0,
#     0,0,0,0,0,0,0,0,0,0,0,0,
#     0,0,0,0,0,0,0,0]
    






# if __name__ == '__main__':

#     # final_gene = [0,0,0,0,0,1,0,1,0,0,0,0,
#     # 0,1,1,0,0,0,0,1,0,1,0,1,
#     # 1,0,0,0,0,1,1,1,0,0,0,0,
#     # 0,1,1,1,0,0,0,0,0,1,0,0,
#     # 0,0,0,1,0,1,0,1,0,1,0,1,
#     # 0,1,1,1,0,1,1,0,0,1,0,0,
#     # 0,1,1,1,0,1,1,1,0,0,0,0,
#     # 0,1,0,1,0,0,0,0,0,0,0,1,
#     # 0,1,1,1,1,1,0,1,1,1,1,1,
#     # 1,1,1,1,1,0,1,1,0,1,1,1,
#     # 0,1,1,1,1,1,1,1]



#     # gene_length = 2**(2*3+1)
#     # gene = np.random.randint(2,size = gene_length)


#     gene_pool_initial = np.loadtxt('gene_pool').astype(int)

#     avr_time = 0
#     loop = 50
#     gene_pool = gene_pool_initial[0:loop]
#     # print(gene_pool)

#     for i in range(loop):
#         stime = time.time()    
#         indiv = individual(gene_pool[i])
#         etime = time.time() - stime
#         avr_time += etime
#     print('循环用时:{}秒'.format(avr_time))

#     stime = time.time()
     
#     #python multiprocessing 并行计算
#     #windows 下只能在main函数里写
#     pool = Pool(5)
#     pool.map(individual, gene_pool)
#     pool.close()
#     pool.join()

#     etime = time.time() - stime
        
#     print('并行用时:{}秒'.format(etime))
 
