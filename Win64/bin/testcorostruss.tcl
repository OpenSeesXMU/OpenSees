# Reference: An improvedperidynamicapproachforquasi-staticelastic deformation andbrittlefractureanalysis
wipe;
model BasicBuilder -ndm 2 -ndf 2

# --------------------------------------------------------------------------------------------------------------
# N O D E S
# --------------------------------------------------------------------------------------------------------------
set L 1.0; set H 0.2;
# global nbound;
# global nL;
# global nH;
set nbound 5;
set nL 50;
set nH 10;

set dx [expr $L/$nL];puts "dx: $dx"
set radij [expr $dx/2.0];set t 1.0;
set horizon [expr 3.015*$dx];
puts "horizon: $horizon"
set volumn [expr $dx*$dx*$t];
puts "volumn: $volumn"
set E 22e9;set v [expr 1.0/3.0];
set pi [expr 2*asin(1.0)];
set Ec [expr 315.0*$E/8.0/$pi/pow(3.015,3)/$t];
# set Ec [expr 6.0*$E/($pi*$t*pow($horizon,3.0)*(1.0-$v))];
set Ed 0.0;


# puts $nL
# puts $nH

set s0 0.05;
# node $NodeTag $XCoord $Ycoord
set nodeID 0
for {set i 1} {$i <= [expr $nL+$nbound]} {incr i} {
	for {set j 1} {$j <= $nH} {incr j} {
	    set nodeID [expr $nodeID+1]
	    node $nodeID [expr ($i-1)*$dx+$radij - $nbound*$dx] [expr ($j-1)*$dx+$radij]
		puts "node $nodeID [expr ($i-1)*$dx+$radij - $nbound*$dx] [expr ($j-1)*$dx+$radij]"
	}
}
puts "total number of node: $nodeID"
# set Cnode 200000;
# node $Cnode $L [expr $H/2.0]
# for {set i [expr $nL]} {$i < [expr $nL+$nbound]} {incr i} {
	# for {set j 1} {$j <= $nH} {incr j} {
		# equalDOF $Cnode [expr $i*$nH+$j] 1 2
		# puts "equalDOF $Cnode [expr $i*($nH+1)+$j] 1 2"
	# }
# }

set trussID 0;
uniaxialMaterial Elastic 1 [expr $Ec/10000] 
for {set inode 1} {$inode < $nodeID} {incr inode} {
	# puts "node:$inode neigbhour list: "
	for {set jnode [expr $inode+1]} {$jnode <= $nodeID} {incr jnode} {
		set nodeICoordX [nodeCoord $inode 1]
		set nodeICoordY [nodeCoord $inode 2]
		set nodeJCoordX [nodeCoord $jnode 1]
		set nodeJCoordY [nodeCoord $jnode 2]
		set LIJ [expr sqrt(($nodeICoordX-$nodeJCoordX)*($nodeICoordX-$nodeJCoordX)+($nodeICoordY-$nodeJCoordY)*($nodeICoordY-$nodeJCoordY))]
		set numNeigh 0;
		
		if {$LIJ < $horizon} {
			set trussID [expr $trussID + 1]
			element corotTruss $trussID  $inode $jnode 1.0 1;
			set numNeigh [expr $numNeigh + 1]
			set neigbList($inode,$numNeigh) $jnode
			# puts -nonewline " $neigbList($inode,$numNeigh)"
		}
		
	}
	# puts -nonewline "\n"
}
puts "total number of element: $trussID"

for {set i 0} {$i < $nbound} {incr i} {
	for {set j 1} {$j <= $nH} {incr j} {
		fix [expr $i*$nH+$j]  1 1
	}
}


set data Out_truss
file mkdir $data

recorder Node -file [format "$data/Node_displacements%i_%i_%i.out" $nbound $nL $nH] -time -nodeRange 1 $nodeID -dof 1 2 disp
recorder Node -file [format "$data/Node_rotations%i_%i_%i.out" $nbound $nL $nH] -time -nodeRange 1 $nodeID -dof 3 disp
recorder Node -file [format "$data/Node_forceReactions%i_%i_%i.out" $nbound $nL $nH] -time -nodeRange 1 $nodeID -dof 1 2 reaction

recorder Node -file [format "$data/Node_momentReactions%i_%i_%i.out" $nbound $nL $nH] -time -nodeRange 1 $nodeID -dof 3 reaction
recorder Element -file [format "$data/Truss_axialForce%i_%i_%i.out" $nbound $nL $nH] -time -eleRange 1 $trussID axialForce
recorder Element -file [format "$data/Truss_deformations%i_%i_%i.out" $nbound $nL $nH] -time -eleRange 1 $trussID deformations
# source display.tcl
printGID truss.msh


logFile "truss_PD.log"
set time_start [clock seconds]
puts "\nAnalysis started  : [clock format $time_start -format %H:%M:%S]"
puts ""
# Loads - Plain Pattern
set P -1000;
pattern Plain 100 Linear {
    
	# for {set i [expr $nL+$nbound-1]} {$i <= [expr $nL+$nbound-1]} {incr i} {
		# for {set j 1} {$j <= $nH} {incr j} {
			# load [expr ($nL+$nbound-1)*$nH+$j] 0 [expr $P/($nbound*$nH)];
		# }
	# }
	for {set j 1} {$j <= $nH} {incr j} {
		load [expr ($nL+$nbound-1)*$nH+$j] 0 [expr $P/$nH];
	}
	# load $nodeID 0 -1000
}    
system BandGeneral
numberer RCM
constraints Transformation
integrator LoadControl 1
test NormDispIncr  1.0e-8 200 2
algorithm Newton
analysis Static

set AnalOk [analyze 1000]

set time_end [clock seconds]
set analysisTime [expr $time_end-$time_start]
puts "Analysis finished : [clock format $time_end -format %H:%M:%S]"
puts "Analysis time     : $analysisTime seconds"
wipe;
