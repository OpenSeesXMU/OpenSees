 wipe ; 
 model basic -ndm 2 -ndf 2
 
 node 1   0.0  0.0
 node 2   1.0  0.0   
 node 3   2.0  0.0
 node 4   3.0  0.0
 node 5   4.0  0.0   
 node 6   5.0  0.0 
 node 7   0.0  1.0
 node 8   1.0  1.0   
 node 9   2.0  1.0
 node 10   3.0  1.0
 node 11   4.0  1.0   
 node 12   5.0  1.0 
 
 fix 1 1 1 
 fix 7 1 1
 
 
 uniaxialMaterial  Elastic 1 20    
 
 element corotTruss2D     1       1     2    1     1
 element corotTruss2D     2       2     3    1     1
 element corotTruss2D     3       3     4    1     1
 element corotTruss2D     4       4     5    1     1
 element corotTruss2D     5       5     6    1     1
 element corotTruss2D     6       7     8    1     1
 element corotTruss2D     7       8     9    1     1
 element corotTruss2D     8       9     10    1     1
 element corotTruss2D     9       10     11    1     1
 element corotTruss2D     10      11     12    1     1
 element corotTruss2D     11      7     2    1     1
 element corotTruss2D     12      8     3    1     1
 element corotTruss2D     13      9     4    1     1
 element corotTruss2D     14      10     5    1     1
 element corotTruss2D     15      11     6    1     1
 element corotTruss2D     16     8    2      1     1
 element corotTruss2D     17     9    3      1     1
 element corotTruss2D     18     10    4      1     1
 element corotTruss2D     19     11    5      1     1
 element corotTruss2D     20     12    6      1     1
 
 
 
 
 recorder Node -file nodeDisp1.out -time -nodeRange  1 12  -dof 1 2 disp
 #recorder Element -file stress11.out  -time  -ele 1  -material stress 
 #recorder Element -file strain11.out  -time  -ele 1  -material strain


 set P -10;
pattern Plain 100 Linear {

	load  12  0  $P;

} 
 
 
 
system BandGeneral
numberer RCM
constraints Transformation
integrator LoadControl 1
test NormDispIncr  1.0e-8 200 2
algorithm Newton
analysis Static

analyze 1000

 