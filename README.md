#my_chamber

###Intro to game programming Project I

The idea for a project was to use Octet's implementation of fluid dynamics to create simple 2D game:
Key points:
  * player moving around a room with movable boxes (by pushing), 
  * boxes block the gas flow
  * aim of the game: to stop the gas from spreading by enclosing it withing area that is bounded by room's walls and boxes. Player is successfull upon "defenging" given percentile of area of whole room. 

Initially I tried to implement to add boxes as sprite class instances adopted from another demo: *example_invaderers*. It turned out that this was not possible as sprite drawing pipline was not meant to use along meshes. To deal with that I have created *mesh_sprite* class, representing a simple mesh consisting of two triangles forming a rectangle. The mesh boundaries were to be mapped to fluid grid and processed.

To simulate fluid-box interaction I've implemented (my_boundary)[https://gist.github.com/witold-gawlowski/3b2db97697c5b3f577355a791d678593] function. Version of the code optimized for one big box worked O.K. but unfortunately any more boxes were fatal for the performance.
Therefore I decided to resign from mesh boxes mapped to fluid grid (which was the expensive part) and represent boxes directly in the grid coordinates.
