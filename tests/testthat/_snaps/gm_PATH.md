# matching correspondence between graph1 and graph2

    gm(A = A, B = B, seeds = seeds, method = "PATH")
    
    Match (10 x 10):
       corr_A corr_B
    1       1      1
    2       2      2
    3       3      3
    4       4      4
    5       5      5
    6       6      6
    7       7      7
    8       8      8
    9       9      9
    10     10     10

---

          [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
     [1,]    1    0    0    0    0    0    0    0    0     0
     [2,]    0    1    0    0    0    0    0    0    0     0
     [3,]    0    0    1    0    0    0    0    0    0     0
     [4,]    0    0    0    1    0    0    0    0    0     0
     [5,]    0    0    0    0    1    0    0    0    0     0
     [6,]    0    0    0    0    0    1    0    0    0     0
     [7,]    0    0    0    0    0    0    1    0    0     0
     [8,]    0    0    0    0    0    0    0    1    0     0
     [9,]    0    0    0    0    0    0    0    0    1     0
    [10,]    0    0    0    0    0    0    0    0    0     1
