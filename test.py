import numpy as np




if __name__ == '__main__':


    for i in range(20):
        adapt = i
        # print(adapt)
        with open('adapt', 'a') as f:
            np.savetxt(f,[adapt])


