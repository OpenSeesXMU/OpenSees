wipe;						            
model basic -ndm 2 -ndf 2;	

node    1            0     0      
node    2            0    0.1  
     
node    3           0.1      0     
node    4            0.1     0.1         




fix  1       1   1   
fix  2       1   1   


 
       


# its okay nDMaterial PlasticityBond2d  11   18500  22  -0.055  3.243e4  66  -0.08  7400  12.54  -0.08
# trying nDMaterial PlasticityBond2d  11   20000  23  -0.055  3.5059e4  69  -0.08  8000  13.11  -0.08
#its perfect nDMaterial PlasticityBond2d  11   20000  23  -0.045  3.5059e4  69  -0.08  8000  13.11  -0.08
#nDMaterial PlasticityBond2d  11   21000  23  -0.04  3.6812e4  69  -0.075  8400  13.11  -0.075
#nDMaterial PlasticityBond2d  11   21000  23  -0.04  3.6812e4  69  -0.075  8400  13.11  -0.075        
nDMaterial TruncatedDP 1 2 2400000 1.4957e+10 1.7677e+10 0.11 2.6614e7 -2.0684e6 1.0e-8	     



element  quad     1  1 2 3 4 1.0 "PlaneStrain" 1 


pattern UniformExcitation    1     1    -accel "Series -factor 1 -filePath acce.txt -dt 0.01"                                 
                                   

 
parameter 11
addToParameter 11 -element 1 material 1 K
addToParameter 11 -element 1 material 2 K
addToParameter 11 -element 1 material 3 K
addToParameter 11 -element 1 material 4 K

updateParameter 11   2e+11 

constraints Transformation;    				                                                              
numberer Plain;					                                                              # 节点编号方法
system BandGeneral;				                                                              # 存储和求解方稿
test RelativeNormUnbalance 0.01 100 2;			                                                              # 收敛性检骿
algorithm Linear;					                                                      # 求解方法采用NEWTON算法
integrator CentralDifference;				                                                      # 荷载控制,定义时间步长 
analysis Transient;                          	                                                            
analyze 3   0.01;










