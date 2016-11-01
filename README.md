#my_chamber

###Intro to game programming Project I

The idea for a project was to use Octet's implementation of fluid dynamics to create simple 2D game:
Key points:
  * player moving around a room with movable boxes (by pushing), 
  * boxes block the gas flow
  * aim of the game: to stop the gas from spreading by enclosing it withing area that is bounded by room's walls and boxes. Player is successfull upon "defenging" given percentile of area of whole room. 

Initially I tried to implement to add boxes as sprite class instances adopted from another demo: *example_invaderers*. It turned out that this was not possible as sprite drawing pipline was not meant to use along meshes. To deal with that I have created 
  * *quad* - abstract mathematical class representing arbitratily transformed rectangle, located in *math* folder.
  * *mesh_sprite* class, representing a simple mesh consisting of two triangles forming a rectangle.
