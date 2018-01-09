# Kinematic_Wave

The Kinematic Wave equation describes one-dimensional steady uniform flow waves moving through 
a stream channel. It is a simplification of the continuity and momentum equations. Kinematic wave 
assumes the pressure term, the convective acceleration, and the local acceleration terms of the momentum
equation to be insignificant, thus only gravity and friction forces are used (Chapra, 1997). 
	Continuity Equation (mass balance): 
█(∂Q/∂x+(∂A_c)/∂t=0 #(1) )
	Momentum Equation:
█(〖g(S〗_o-S_f)=0#(2) )
The Manning’s Equation can now be used by substituting S_o=S_f and R=A_c/P. 
█(Q=1/n  (A_c^(3/5))/P^(2/3)  √(S_o )#(3) )
Then the Manning’s Equation can be solved for:
█(A_c=αQ^β#(4) )
Where β=3/5 and α=〖((nP^(2/3))/√( S_o ))〗^(3/5)
Equation (1) can now be substituted with equation (4) differentiated with respect to time:
█(∂Q/∂x+αβQ^(β-1)  ∂Q/∂t=0 #(5) )
Using a numerical solution by substituting forward time/ backward space:
█((Q_(i,t)-Q_(i-1,t))/∆x+ αβ〖〖(Q〗_(i,t))〗^(β-1)  (Q_(i,t+1)-Q_(i,t))/∆t#(6) )
Where Q_(i,t) is the outflow from the preceding time step, Q_(i-1,t)  is the upstream flow, Q_(i,t+1) is the outflow for this time step.
Then solving for the unknown outflow for this time step:
█(Q_(i,t+1)=1/(αβQ_(i,t)^(β-1) )*(Q_(i-1,t)-Q_(i,t))/∆x ∆t+Q_(i,t)#(7) )
This equation assumes a minimum initial flow in all segments (initial flow ≠ 0)



 
Alpha, α , will change with the shape of the channel
 


 
Channel Shape	Wetted Perimeter	Alpha
Rectangular	2d+w	(n/√(S_o ) (2d+w)^(2/3) )^(3/5)
Trapezoidal	w+2d√(1+z^2 )	(n/√(S_o ) (w+2d√(1+z^2 ))^(2/3) )^(3/5)
Triangular	2d√(1+z^2 )	(n/√(S_o ) (2d√(1+z^2 ))^(2/3) )^(3/5)
 
 
Unchanging parameters for stream segment:
Manning’s roughness	n
Friction slope		S_f
Length segment		L	(m)
Bottom Width		w	(m)
Beta			β	(3/5)
Z-Slope			z	(2)
Changing variables:
	Depth		d	(m)
	Volume		V	(m3)
	Initial Velocity	v_i	(m/s)
	Discharge	 Q	(m3/s)
	Alpha		α

References: 

Chapra, S. C. (1997). Surface Water-Quality Modeling: The McGraw-Hill Companies, Inc. .
 
