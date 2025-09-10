# Grid Generator
This will bet he repository for Aneesh's work on deveoping a grid generator, which will make 3D curvilinear meshes sufficiently smooth to use in numerical simulations.

## Features
Features ofthe grid generator will be:
- Metric for smoothness
- Optimal smoothing algorithm
- User interface
- Ability to refine the grid in a region ( internal points, along edges ... )

## Papers
Included are some papers about grid generation, intermediat grids, and other metrics and methods for grid generation.
- (High Order curvilinear mesh generation from third party meshes)[papers/Curvilinear_Mesh-main.pdf]
  - Included Jacobian metrics for goodness of mesh, and bins them based on valeu as metric
- (On the folding of numerically generated grids)[papers/Folding_Numerically_Generated_Grids1988.pdf]
  - Algorithms to minimize or cancel grid folding.
- (Mathematical Aspects of variational grid generation II)[papers/MathGridGenII.pdf]
  - Original Castillo paper on grid generation using variational grid, using a reference grid
- (Solution adaptive direct variational grids for fluid flow calculations)[papers/Solution_adaptive_direct_variational_grids_for_flu.pdf]
  - Adaptive methods for different grid resolution and marriage of solutions
- (A generalized length startegy for direct optomization in planar grid generation)[papers/a-generalized-length-strategy-for-direct-optimization-in-planar-grid-generation.pdf]    
