 wipe ; 
 model basic -ndm 2 -ndf 2
 node 1   0.0  0.0
 node 2   1.0  0.0   
 node 3   2.0  0.0
 node 4   3.0  0.0
 node 5   4.0  0.0   
 node 6   5.0  0.0 
 node 1   0.0  1.0
 node 2   1.0  1.0   
 node 3   2.0  1.0
 node 4   3.0  1.0
 node 5   4.0  1.0   
 node 6   5.0  1.0 
 
 uniaxialMaterial testSteel01     1   2e1    3.0e-2   0.2   
 
 element corotationalTruss2D     1       1     2    1     1
 element corotationalTruss2D     2       2     3    1     1
 element corotationalTruss2D     3       3     4    1     1
 element corotationalTruss2D     4       4     5    1     1
 element corotationalTruss2D     5       5     6    1     1
 element corotationalTruss2D     7       7     8    1     1
 element corotationalTruss2D     8       8     9    1     1
 element corotationalTruss2D     9       9     10    1     1
 element corotationalTruss2D     10      10     11    1     1
 element corotationalTruss2D     11      11     12    1     1
 element corotationalTruss2D     12      7     2    1     1
 element corotationalTruss2D     13      8     3    1     1
 element corotationalTruss2D     14      9     4    1     1
 element corotationalTruss2D     15      10     5    1     1
 element corotationalTruss2D     16      11     6    1     1
 
 
 recorder Node -file nodeDisp.out -time -nodeRanage 1 12  -dof 1 2 disp
 #recorder Element -file stress11.out  -time  -ele 1  -material stress 
 #recorder Element -file strain11.out  -time  -ele 1  -material strain


 set P -1000;
pattern Plain 100 Linear {

	load  12  0 $P;

} 
 
 
 
 constraints Transformation
 numberer RCM
 test NormDispIncr 1.E-8 25  2
 algorithm Newton
 system BandSPD
 integrator Newmark 0.55 0.275625 
 analysis Transient
 analyze 1000  0.01

 