# perco of same sizes

    graph_match_percolation(A = A, B = B, seeds = seeds)
    
    Match (10 x 10):
      corr_A corr_B
    1      1      1
    2      3      3
    3      4      8
    4      5      5
    5      6      6
    6      7      7
    7      9      9
    8     10      2

---

    WAoAAAACAAQAAwACAwAAAQMTAAAAAgAAAA0AAAAIAAAAAQAAAAMAAAAHAAAAAgAAAAYAAAAF
    AAAABAAAAAgAAAMTAAAAAgAAAA4AAAADP/AAAAAAAABAFAAAAAAAAEAIAAAAAAAAAAAADgAA
    AAM/8AAAAAAAAEAUAAAAAAAAQAgAAAAAAAAAAAQCAAAAAQAEAAkAAAAFbmFtZXMAAAAQAAAA
    AgAEAAkAAAABQQAEAAkAAAABQgAABAIAAAABAAQACQAAAAVjbGFzcwAAABAAAAABAAQACQAA
    AApkYXRhLmZyYW1lAAAEAgAAAAEABAAJAAAACXJvdy5uYW1lcwAAAA0AAAACgAAAAP////0A
    AAD+AAAEAgAAAf8AAAAQAAAAAgAEAAkAAAALbWF0Y2hfb3JkZXIABAAJAAAABXNlZWRzAAAE
    AgAAAAEABAAJAAAABGNvcnIAAAMTAAAAAgAAAA4AAAAIP/AAAAAAAABACAAAAAAAAEAQAAAA
    AAAAQBQAAAAAAABAGAAAAAAAAEAcAAAAAAAAQCIAAAAAAABAJAAAAAAAAAAAAA4AAAAIP/AA
    AAAAAABACAAAAAAAAEAgAAAAAAAAQBQAAAAAAABAGAAAAAAAAEAcAAAAAAAAQCIAAAAAAABA
    AAAAAAAAAAAABAIAAAH/AAAAEAAAAAIABAAJAAAABmNvcnJfQQAEAAkAAAAGY29ycl9CAAAE
    AgAAA/8AAAAQAAAACAAEAAkAAAABMQAEAAkAAAABMgAEAAkAAAABMwAEAAkAAAABNAAEAAkA
    AAABNQAEAAkAAAABNgAEAAkAAAABNwAEAAkAAAABOAAABAIAAAL/AAAAEAAAAAEABAAJAAAA
    CmRhdGEuZnJhbWUAAAD+AAAEAgAAAAEABAAJAAAABm5ub2RlcwAAAA0AAAACAAAACgAAAAoA
    AAQCAAAAAQAEAAkAAAAEY2FsbAAAAAYAAAABAAQACQAAABdncmFwaF9tYXRjaF9wZXJjb2xh
    dGlvbgAABAIAAAABAAQACQAAAAFBAAAI/wAABAIAAAABAAQACQAAAAFCAAAJ/wAABAIAAAAB
    AAQACQAAAAVzZWVkcwAACv8AAAD+AAAEAgAAAv8AAAIQAAAAAQAEAAkAAAAKZ3JhcGhNYXRj
    aAAABAIAAAABAAQACQAAAAdwYWNrYWdlAAAAEAAAAAEABAAJAAAAC2lHcmFwaE1hdGNoAAAA
    /gAAAP4=

# perco w. similarity score

    graph_match_percolation(A = A, B = B, seeds = seeds, similarity = sim)
    
    Match (10 x 10):
      corr_A corr_B
    1      1      1
    2      2      6
    3      3      3
    4      5      5
    5      8      9
    6      9      7
    7     10      4

---

    WAoAAAACAAQAAwACAwAAAQMTAAAAAgAAAA0AAAAHAAAAAQAAAAcAAAADAAAAAgAAAAQAAAAF
    AAAABgAAAxMAAAACAAAADgAAAAM/8AAAAAAAAEAUAAAAAAAAQAgAAAAAAAAAAAAOAAAAAz/w
    AAAAAAAAQBQAAAAAAABACAAAAAAAAAAABAIAAAABAAQACQAAAAVuYW1lcwAAABAAAAACAAQA
    CQAAAAFBAAQACQAAAAFCAAAEAgAAAAEABAAJAAAABWNsYXNzAAAAEAAAAAEABAAJAAAACmRh
    dGEuZnJhbWUAAAQCAAAAAQAEAAkAAAAJcm93Lm5hbWVzAAAADQAAAAKAAAAA/////QAAAP4A
    AAQCAAAB/wAAABAAAAACAAQACQAAAAttYXRjaF9vcmRlcgAEAAkAAAAFc2VlZHMAAAQCAAAA
    AQAEAAkAAAAEY29ycgAAAxMAAAACAAAADgAAAAc/8AAAAAAAAEAAAAAAAAAAQAgAAAAAAABA
    FAAAAAAAAEAgAAAAAAAAQCIAAAAAAABAJAAAAAAAAAAAAA4AAAAHP/AAAAAAAABAGAAAAAAA
    AEAIAAAAAAAAQBQAAAAAAABAIgAAAAAAAEAcAAAAAAAAQBAAAAAAAAAAAAQCAAAB/wAAABAA
    AAACAAQACQAAAAZjb3JyX0EABAAJAAAABmNvcnJfQgAABAIAAAP/AAAAEAAAAAcABAAJAAAA
    ATEABAAJAAAAATIABAAJAAAAATMABAAJAAAAATQABAAJAAAAATUABAAJAAAAATYABAAJAAAA
    ATcAAAQCAAAC/wAAABAAAAABAAQACQAAAApkYXRhLmZyYW1lAAAA/gAABAIAAAABAAQACQAA
    AAZubm9kZXMAAAANAAAAAgAAAAoAAAAKAAAEAgAAAAEABAAJAAAABGNhbGwAAAAGAAAAAQAE
    AAkAAAAXZ3JhcGhfbWF0Y2hfcGVyY29sYXRpb24AAAQCAAAAAQAEAAkAAAABQQAACP8AAAQC
    AAAAAQAEAAkAAAABQgAACf8AAAQCAAAAAQAEAAkAAAAFc2VlZHMAAAr/AAAEAgAAAAEABAAJ
    AAAACnNpbWlsYXJpdHkAAAABAAQACQAAAANzaW0AAAD+AAAEAgAAAv8AAAIQAAAAAQAEAAkA
    AAAKZ3JhcGhNYXRjaAAABAIAAAABAAQACQAAAAdwYWNrYWdlAAAAEAAAAAEABAAJAAAAC2lH
    cmFwaE1hdGNoAAAA/gAAAP4=

# percolation without seeds

    graph_match_percolation(A = A, B = B, seeds = NULL, similarity = sim)
    
    Match (10 x 10):
      corr_A corr_B
    1      1      2
    2      2      6
    3      3      9
    4      4     10
    5      5      7
    6      6      3
    7      7      4
    8      9      5

---

    WAoAAAACAAQAAwACAwAAAQMTAAAAAgAAAA0AAAAIAAAABgAAAAcAAAAEAAAACAAAAAEAAAAC
    AAAAAwAAAAUAAAMTAAAAAgAAAA4AAAAAAAAADgAAAAAAAAQCAAAAAQAEAAkAAAAFbmFtZXMA
    AAAQAAAAAgAEAAkAAAABQQAEAAkAAAABQgAABAIAAAABAAQACQAAAAVjbGFzcwAAABAAAAAB
    AAQACQAAAApkYXRhLmZyYW1lAAAEAgAAAAEABAAJAAAACXJvdy5uYW1lcwAAAA0AAAAAAAAA
    /gAABAIAAAH/AAAAEAAAAAIABAAJAAAAC21hdGNoX29yZGVyAAQACQAAAAVzZWVkcwAABAIA
    AAABAAQACQAAAARjb3JyAAADEwAAAAIAAAANAAAACAAAAAEAAAACAAAAAwAAAAQAAAAFAAAA
    BgAAAAcAAAAJAAAADQAAAAgAAAACAAAABgAAAAkAAAAKAAAABwAAAAMAAAAEAAAABQAABAIA
    AAH/AAAAEAAAAAIABAAJAAAABmNvcnJfQQAEAAkAAAAGY29ycl9CAAAEAgAAA/8AAAAQAAAA
    CAAEAAkAAAABMQAEAAkAAAABMgAEAAkAAAABMwAEAAkAAAABNAAEAAkAAAABNQAEAAkAAAAB
    NgAEAAkAAAABNwAEAAkAAAABOAAABAIAAAL/AAAAEAAAAAEABAAJAAAACmRhdGEuZnJhbWUA
    AAD+AAAEAgAAAAEABAAJAAAABm5ub2RlcwAAAA0AAAACAAAACgAAAAoAAAQCAAAAAQAEAAkA
    AAAEY2FsbAAAAAYAAAABAAQACQAAABdncmFwaF9tYXRjaF9wZXJjb2xhdGlvbgAABAIAAAAB
    AAQACQAAAAFBAAAI/wAABAIAAAABAAQACQAAAAFCAAAJ/wAABAIAAAABAAQACQAAAAVzZWVk
    cwAAAP4AAAQCAAAAAQAEAAkAAAAKc2ltaWxhcml0eQAAAAEABAAJAAAAA3NpbQAAAP4AAAQC
    AAAC/wAAAhAAAAABAAQACQAAAApncmFwaE1hdGNoAAAEAgAAAAEABAAJAAAAB3BhY2thZ2UA
    AAAQAAAAAQAEAAkAAAALaUdyYXBoTWF0Y2gAAAD+AAAA/g==

# perco w. directed graphs

    graph_match_percolation(A = A, B = B, seeds = seeds, similarity = sim)
    
    Match (10 x 10):
       corr_A corr_B
    1       1      1
    2       2      8
    3       3      3
    4       4      4
    5       5      5
    6       6      6
    7       7      7
    8       8      9
    9       9      2
    10     10     10

---

    WAoAAAACAAQAAwACAwAAAQMTAAAAAgAAAA0AAAAKAAAAAQAAAAoAAAADAAAABAAAAAIAAAAF
    AAAABgAAAAcAAAAIAAAACQAAAxMAAAACAAAADgAAAAM/8AAAAAAAAEAUAAAAAAAAQAgAAAAA
    AAAAAAAOAAAAAz/wAAAAAAAAQBQAAAAAAABACAAAAAAAAAAABAIAAAABAAQACQAAAAVuYW1l
    cwAAABAAAAACAAQACQAAAAFBAAQACQAAAAFCAAAEAgAAAAEABAAJAAAABWNsYXNzAAAAEAAA
    AAEABAAJAAAACmRhdGEuZnJhbWUAAAQCAAAAAQAEAAkAAAAJcm93Lm5hbWVzAAAADQAAAAKA
    AAAA/////QAAAP4AAAQCAAAB/wAAABAAAAACAAQACQAAAAttYXRjaF9vcmRlcgAEAAkAAAAF
    c2VlZHMAAAQCAAAAAQAEAAkAAAAEY29ycgAAAxMAAAACAAAADgAAAAo/8AAAAAAAAEAAAAAA
    AAAAQAgAAAAAAABAEAAAAAAAAEAUAAAAAAAAQBgAAAAAAABAHAAAAAAAAEAgAAAAAAAAQCIA
    AAAAAABAJAAAAAAAAAAAAA4AAAAKP/AAAAAAAABAIAAAAAAAAEAIAAAAAAAAQBAAAAAAAABA
    FAAAAAAAAEAYAAAAAAAAQBwAAAAAAABAIgAAAAAAAEAAAAAAAAAAQCQAAAAAAAAAAAQCAAAB
    /wAAABAAAAACAAQACQAAAAZjb3JyX0EABAAJAAAABmNvcnJfQgAABAIAAAP/AAAAEAAAAAoA
    BAAJAAAAATEABAAJAAAAATIABAAJAAAAATMABAAJAAAAATQABAAJAAAAATUABAAJAAAAATYA
    BAAJAAAAATcABAAJAAAAATgABAAJAAAAATkABAAJAAAAAjEwAAAEAgAAAv8AAAAQAAAAAQAE
    AAkAAAAKZGF0YS5mcmFtZQAAAP4AAAQCAAAAAQAEAAkAAAAGbm5vZGVzAAAADQAAAAIAAAAK
    AAAACgAABAIAAAABAAQACQAAAARjYWxsAAAABgAAAAEABAAJAAAAF2dyYXBoX21hdGNoX3Bl
    cmNvbGF0aW9uAAAEAgAAAAEABAAJAAAAAUEAAAj/AAAEAgAAAAEABAAJAAAAAUIAAAn/AAAE
    AgAAAAEABAAJAAAABXNlZWRzAAAK/wAABAIAAAABAAQACQAAAApzaW1pbGFyaXR5AAAAAQAE
    AAkAAAADc2ltAAAA/gAABAIAAAL/AAACEAAAAAEABAAJAAAACmdyYXBoTWF0Y2gAAAQCAAAA
    AQAEAAkAAAAHcGFja2FnZQAAABAAAAABAAQACQAAAAtpR3JhcGhNYXRjaAAAAP4AAAD+

# percolation multi-layer

    graph_match_percolation(A = A, B = B, seeds = seeds)
    
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

    WAoAAAACAAQAAwACAwAAAQMTAAAAAgAAAA0AAAAKAAAAAQAAAAIAAAADAAAACQAAAAUAAAAH
    AAAABAAAAAoAAAAGAAAACAAAAxMAAAACAAAADQAAAAMAAAABAAAAAgAAAAMAAAANAAAAAwAA
    AAEAAAACAAAAAwAABAIAAAABAAQACQAAAAVuYW1lcwAAABAAAAACAAQACQAAAAFBAAQACQAA
    AAFCAAAEAgAAAAEABAAJAAAABWNsYXNzAAAAEAAAAAEABAAJAAAACmRhdGEuZnJhbWUAAAQC
    AAAAAQAEAAkAAAAJcm93Lm5hbWVzAAAADQAAAAKAAAAA/////QAAAP4AAAQCAAAB/wAAABAA
    AAACAAQACQAAAAttYXRjaF9vcmRlcgAEAAkAAAAFc2VlZHMAAAQCAAAAAQAEAAkAAAAEY29y
    cgAAAxMAAAACAAAADQAAAAoAAAABAAAAAgAAAAMAAAAEAAAABQAAAAYAAAAHAAAACAAAAAkA
    AAAKAAAADQAAAAoAAAABAAAAAgAAAAMAAAAEAAAABQAAAAYAAAAHAAAACAAAAAkAAAAKAAAE
    AgAAAf8AAAAQAAAAAgAEAAkAAAAGY29ycl9BAAQACQAAAAZjb3JyX0IAAAQCAAAD/wAAABAA
    AAAKAAQACQAAAAExAAQACQAAAAEyAAQACQAAAAEzAAQACQAAAAE0AAQACQAAAAE1AAQACQAA
    AAE2AAQACQAAAAE3AAQACQAAAAE4AAQACQAAAAE5AAQACQAAAAIxMAAABAIAAAL/AAAAEAAA
    AAEABAAJAAAACmRhdGEuZnJhbWUAAAD+AAAEAgAAAAEABAAJAAAABm5ub2RlcwAAAA0AAAAC
    AAAACgAAAAoAAAQCAAAAAQAEAAkAAAAEY2FsbAAAAAYAAAABAAQACQAAABdncmFwaF9tYXRj
    aF9wZXJjb2xhdGlvbgAABAIAAAABAAQACQAAAAFBAAAI/wAABAIAAAABAAQACQAAAAFCAAAJ
    /wAABAIAAAABAAQACQAAAAVzZWVkcwAACv8AAAD+AAAEAgAAAv8AAAIQAAAAAQAEAAkAAAAK
    Z3JhcGhNYXRjaAAABAIAAAABAAQACQAAAAdwYWNrYWdlAAAAEAAAAAEABAAJAAAAC2lHcmFw
    aE1hdGNoAAAA/gAAAP4=
