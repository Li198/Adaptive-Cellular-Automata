


class config:
    dict = {#individual
            "width":49,   
            "height":49,   
            "neighbourhood":3,            # __neighbourhood = get("neighbourhood") #规则同上一层的左右三个细胞有关
            "individual_test_step":100,    #individual 中通过多少步来测试 adaptability
            
            "neighbour_radius":1,
            "group_width":10,             #generation 中CGA的长宽
            "group_height":10,            # "parent_number":200,          
            #generation 中一代有多少父母个体
            # "child_number":200,  
            "generation_num":10,           #generation 跑多少代
            }

    def get(self,key):
        return self.dict[key]

# if __name__ == '__main__':
#     width = config().get("neighbourhood")

#     print(width)
