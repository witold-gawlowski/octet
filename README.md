#my_chamber

###Intro to game programming Project I

> This project is a modification of Andy Thomasson's implementation of well known Jos Stams's paper about fluid mechanics.
> In the original version the fundamental gas propagation is implemented. My work expands on interaction with boxes.
> Specifically it adds:
> * Fluid Drag Begind the object
> * Pushing fluid in from of a moving object
> * Box random Generation.
> * Player character, character movement, health decreas upon gas interaction
> * Pushing boxes around with 
> * Blocking fluid with box volumes.
> * Visual feedback on "leaks": unblocked paths between the player square
> * Simple GUI
> * Score calculation

## 
I have tracked the progress of the project in "outline" file available [here](https://github.com/witold-gawlowski/octet/blob/Intro_to_game_programming_1/octet/src/projects/my_chamber/Concept/Outline.txt).

## Player-box interaction
Box grid is an array that stores the index of a box occupying a given position. When players trys to walk over it, the function move_box with boxes index is called triggering 
all necessary changes. 

## Visual feedback for leaks.
I determine if the player is separated from the fluid source, as well as the separated area, using depth first seach algorithm (DFS). With DFS implemented it was easy to visualize
a path that DFS travels when it reches the source. This is sometimes helpfull when its hard to spot the place where the created bareer is leaking. 

## Fluid-box interaction
To simulate the interaction of dynamic elements with gas I have implemented several interaction mechanisms:
1. Fluid rebound: I have adopted the rebound from chamber's barriers to the dynamic barriers of boxes.
2. Fluid drag behind the boxes: I add manually some values to the velocity fields behind player.
3. Gas accumulation in front of a moving object: I "push" the values of the fluid intensity in front of a moving cube.

[Here](https://github.com/witold-gawlowski/octet/blob/Intro_to_game_programming_1/octet/src/projects/my_chamber/my_chamber.h) is the source file with the fluid-box interaction implementation. I did my best to comment the core parts as clrearly as possible. 

## Sprite class
To ease the work with boxes I have also implemented a Sprite class that is representing a drawable, square mesh, consisting of two triangles as well as a correspoing representaion of unrelying mathematical object defining the geometry. Boxes visible in the game are drawn using that class. It was implemented to fit the structure of other drawable calsses in Octet. Code of the mathematical class implemented is available [here](https://github.com/witold-gawlowski/octet/blob/Intro_to_game_programming_1/octet/src/projects/my_chamber/quad.h). 



## Game mechanics

I have arranged implemented fluid interaction into a simple game where player is to stop the biggest possible area from being polluted while staying alive
(player looses health according to pollution intensity at the given point). Score is calculated as

> (area separated from fluid source)/(polution in separated area + 1).

I add 1 to the denominator to prevent division by 0 and also, in case of no pollution the score is just area separated from fluid.
Player is assigned zero score if there is a path from fluid source to the player using non-diagonal edges.

## Summary
The project started to went on really quick when I have solved some performance problems caused by using vector overlay over data (instead of accessing the data directly).
I have greatly enjoyed adopinng the simulation for interaction, Jos Stams paper as well as the design of the little game's mechanics. Overall it was a great project to work on!

