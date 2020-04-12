 wipe ; 
 model basic -ndm 2 -ndf 2
 node 1   0.0  0.0
 node 2   10.0  0.0   -mass 10000.0  10000.0
 fix 1 1 1 
 fix 2 0 1
 uniaxialMaterial testSteel01     1   2e1    3.0e-2   0.2   
 uniaxialMaterial Elastic   2   20.0
 element truss     1      1     2    1     1
 element truss     2      1     2    1     2
 recorder Node -file node2.out -time -node 2 -dof 1 2 disp
 recorder Element -file stress1.out  -time  -ele 1  -material stress 
 recorder Element -file strain1.out  -time  -ele 1  -material strain
 recorder Element -file stress2.out  -time  -ele 2  -material stress 
 recorder Element -file strain2.out  -time  -ele 2  -material strain
 set tabas "Path -filePath tabas.txt -dt 0.02 -factor 4"
 pattern UniformExcitation  1   1  -accel      $tabas
 constraints Transformation
 numberer RCM
 test NormDispIncr 1.E-8 25  2
 algorithm Newton
 system BandSPD
 integrator Newmark 0.55 0.275625 
 analysis Transient
 analyze 1000  0.01
