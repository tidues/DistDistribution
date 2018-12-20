import numpy as np

class array1d:
    # init an indexed array
    def __init__(self, keyLst, valLst):
        # build
        self.arr = np.array(valLst)
        self.keyLst = keyLst
        # check shape errors
        if len(keyLst) != len(valLst):
            print("keyLst shape is incorrect")
            return -1

        # build dict
        self.keyMap = {}
        for idx,key in enumerate(keyLst):
            self.keyMap[key] = idx

    # define slicing by index
    def __getitem__(self, key):
        return self.arr[self.keyMap[key]]

class array2d:
    def __init__(self, keysLst, arrLst):


if __name__ == '__main__':
    
    lst0 = [0,1,2]
    lst1 = [3,4,5]
    keys = ['c','d','e']
    arr0 = array(keys, lst)
    arr1 = array(keys, lst)
    print(arr['c'])
    print(arr['d'])
    keys1 = ['a','b']
    arr2d = array(keys1, [arr0, arr1])



    
