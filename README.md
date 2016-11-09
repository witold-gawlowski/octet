# my_bridge: Tools and Middleware Project I

Final project code Features:
* Collada mesh of the canyon created by Martin Skarregard.
* Bullet physics bridge:
	* bt

The core goal of the project was to create a bullet-phyiscs demo. Our lecturer - Andy Thomason suggested building a rope bridge. I liked the idea as it fitted the concept one of my projects I wanted to realize: an experience of a canyon: I have found this picture somewhere on the web and wanted to make a game where you could explore it. 

## Final project Features:
* Displaying Collada mesh of the canyon created by Martin Skarregard.
* Bullet physics bridge.
* Octet FPS controller
* Fullscreen, no title bar mode
* bpy (blender API for python) [script](https://github.com/witold-gawlowski/octet/blob/Tools_and_Middleware_1/octet/assets/projects/my_bridge/scripts/box_tracker.py) that reads the scene and writes the bridge and world colliders elements to the csv file
* reading world colliders and bridge elements from csv file in Octet code

## Goals and inspirations 

The main goal of the project was to create a bullet-phyiscs demo. We were given the idea of builing a rope bridge. I liked it as it fitted the frames of of one of the projects I wanted to realize: an experience of a canyon. I have come a cross the idea for a canyon projecet when due to following graphics that I've found on the internet and imagined to be nice idea to make a game where you could explore it. 

![inspirtation](https://github.com/witold-gawlowski/octet/blob/Tools_and_Middleware_1/octet/src/projects/my_bridge/Concept/14446003_1212966025441859_1981789854985566161_n.jpg) 

For that I needed a canyon mesh. Initially I wanted to create it myself and learn blender as a side-effect but decided to ask on the facebook group first. Martin Skarregaard has answered me and we started to collaborate. I have created a story-board for him: 

![](https://github.com/witold-gawlowski/octet/blob/Tools_and_Middleware_1/octet/src/projects/my_bridge/Concept/canyon%20mesh%20specification.jpg)

I was happy with Martins work although many things didn't match my specification. Also some things that I thought were obvious were not included.
 
I have created much less content for this project than I have initially planned to. Informal backlog featuring ideas and technical problem encountered during the project can be found in [issue_track.txt file](https://github.com/witold-gawlowski/octet/blob/Tools_and_Middleware_1/octet/src/projects/my_bridge/Concept/issue_tracker.txt). There are two reasons for this: 


i got lost in the second assignement

After finishing my second assigment i had 2 days to make more interesting rope-brigdge. I have reworked the wole code architecture, made it object oriented (bridge was to be consisting of multiple instances of segment class)

Although the new code architecture worked properly, the physics-engine problems were very severe, the parts were bouncing and wobbling for the reason i have not yet found out. I have started to build same bridge directly in bullet3's examples code and it has been working so far. This gives me a starting point for further investigation. 

Problems I encoutered:
* Bullet library is poorly documented and its not easy to find information on the internet too. One has to basically undestand the implementation of api to be able to use it proficiently.
* Using physics engines is difficult(with a slight chance that I have been just doing something really wrong).


What I have learned:
* When you specify a task for someone you need to point out every detail that is significant. Ex. If you want proportions on the model to be 3:4 you need it point it out, not just draw it 3:4.

The video with the outcome of our work could be found [here][https://www.youtube.com/watch?v=uHrpnjAkioI]. I like the atmosphere result of the project although I am unstisfied with the amout of features I was able to implement. The elements missing the most for me are:
* fmod sound effects
* skybox with stars
* colliders for the walls of the caves
* 
=======
I was happy with Martins work, although many things didn't match my specification. Also some things that I thought were obvious were not included.
The Canyon he scuplted is the one visible on the video linked later in this readme. 

## Post-mortem
 
I have created much less content for this project than I have initially planned to. Informal backlog for feature ideas and some technical problem encountered during the project can be found in [issue_track.txt file](https://github.com/witold-gawlowski/octet/blob/Tools_and_Middleware_1/octet/src/projects/my_bridge/Concept/issue_tracker.txt). There are several reasons for that. I have prioritized  the parrallel assignment for Intro to Game Programming over this one. Also, my goal for further development for this project was to create a "real rope bridge", but some technical problems with bullet arose. 

Initially I have reworked the code architecture, made it object oriented (bridge was to be consisting of multiple instances of segment class) that is visible on the tip of the branch. 
Although the new code architecture worked properly, the physics-engine problems were very severe, the parts were bouncing and wobbling for the reason I have not yet found out. I have started to build simmilar bridge directly in bullet3's example code and it has been working so far. This gives me a starting point for further investigation. 

## Some problems I have encountered:
* Bullet library is poorly documented and its not easy to find information on the internet too. One has to basically undestand the implementation of api to be able to use it proficiently.
* Using physics engines is difficult (with a slight chance that I have been just doing something really wrong).

## What I have learned:
* Python-blender scripting is extremly effective for making tools. It is also convenient. E.g.: python API visible when hovering over any blender GUI element. 
* Writing code comments about what i have learned about library i am using is important. Its not easy to remember the functional definition for some variable, that you have come up after long time of reversed-engineering.
* During my attempts to build more complex rope bridge I mastered 90% of btGeneric6DOFConstraint API's which I feel very proud of. 
* When you specify a task you need to point out every detail that is significant. Ex. If you want proportions on the model to be 3:4 you need it point it out, not just draw it 3:4 and this is also a lot of work so unless you have a reason its better to do it yourself. 

## The video with the outcome of our work could be found [here][https://www.youtube.com/watch?v=uHrpnjAkioI]. I like the atmosphere result of the project although I am unstisfied with the amount of features I was able to implement. The elements missing the most for me are:
* fmod sound effects (the sounds in youtube vidoe were added throught video editing)
* skybox with stars
* colliders for the walls of the caves
* custom shader for canyon (faking ambient occlusion for the bottom of the canyon)
* lightmapping 
* vertex shader for canyon to make it look more rough
* more complex bridge physics

Hopefully I will find some time to add some of those for my portfolio. 

The commit, that I would like to be marked upon is [here](https://github.com/witold-gawlowski/octet/tree/e4f9d0dff4dcd4b0e8cb072813f0d4e0aed31ad5).
The tip of the branch, althought is not working, includes a lot of work: refactored, object oriented code among other things. The thing that is not working is physics.

