### normals 
Rendering is weird with some black and white triangles
* check the order of the vertices when computing the normal (i.e. v1-v2 and v1-v3)
* check cross is used correctly --> it returns a value, it's not a modifier

### no render
Nothing is shown (especially in the early stage)
* check the data collection is ok
* check they subtract -1 to the indices
* check they didn't declare some local variable where to store the data inside load()
* check the path is correct

### index rendering
* check that they are giving the good number of elements