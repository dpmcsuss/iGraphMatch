# graphMatch class functions

    Call: gm(A = g1, B = g2, method = "indefinite", start = "bari")
    
    # Matches: 10
    # True Matches:  0, # Seeds:  0, # Vertices:  10, 200
                    
      common_edges 7
     missing_edges 2
       extra_edges 0
             fnorm 2

# as.graphMatch works

    Code
      as.graphMatch(data.frame(a = 1:5, b = 1:5))
    Output
      asMethod(from = object)
      
      Match (5 x 5):
        corr_A corr_B
      1      1      1
      2      2      2
      3      3      3
      4      4      4
      5      5      5

---

    Code
      as.graphMatch(data.frame(corr_A = 1:5, corr_B = 1:5))
    Output
      asMethod(from = object)
      
      Match (5 x 5):
        corr_A corr_B
      1      1      1
      2      2      2
      3      3      3
      4      4      4
      5      5      5

---

    Code
      as.graphMatch(1:10)
    Output
      asMethod(from = object)
      
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

