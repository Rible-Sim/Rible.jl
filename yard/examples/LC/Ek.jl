qs = sol_freevibra001.qs
len = length(qs)
r3 = [[qs[i][5],qs[i][6]] for i = 1:len]
r4 = [[qs[i][7],qs[i][8]] for i = 1:len]
r5 = [[r4[i][2]-r3[i][2],r3[i][1]-r4[i][1]] for i = 1:len] + r4;
r1 = [[0;0] for i = 1:len ]
r2 = [[-0.1;0.1] for i = 1:len]
r6 = [[0.0;0.2] for i = 1:len ]
ul = 0.01norm(r3[1]-r2[1])
ll = 0.01norm(r3[1]-r1[1])
K = (norm.(r3 - r1).-ll).^2+(norm.(r3 - r2).-ul).^2+(norm.(r5 - r2).-ul).^2+(norm.(r5 - r6).-ll).^2
