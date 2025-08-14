Bouncing Ball
=============
   
.. raw:: html

This example graphically illustrates the loss/gain of potential energy over time due to damping effects as the simulation is run forward/backward in time. If the ball were to ever reach a steady state condition where it would fully come to rest, then it would be theoretically impossible to reverse the flow of time and determine the initial starting position of the ball if all that one knew was the terminal (steady state) condition of the ball. However, given that the proposed methodology keeps track of the energy incrementally discarded over the full duration of a particular forward-in-time analysis, it is nonetheless possible within this framework to exactly reverse this analysis to exactly recover the initial conditions. It is nonetheless important to emphasize that without this ancillary information regarding the continuous dissipation of energy (information), there are many possible prior trajectories that the state of the system might have followed to achieve the terminal steady state conditions. Hence, it is not appropriate to reverse the flow of time directly given starting terminal conditions.

	 <iframe src="../_static/bouncing_ball/index.html" style="width:100%; height:100%; aspect-ratio: 10/6;" scrolling="no" frameborder="0"></iframe>
