An Augmented Reality application
===========================================

----------------------
Building instructions
----------------------

Required tools:
* C/C++ compiler (gcc >= 4.6 or visual studio or clang)
* cmake
* [optional] Doxygen and graphviz to generate the htlm documentation


###  Dependencies

The project depends on:

- OpenGL


If at any stage you get an error about missing ``Xi`` and ``Xmu`` library, you need to install them. In linux you can do

```shell
sudo apt-get install libxi-dev libxmu-dev.
```

###  Setting up and building:

OpenGL is normally already installed on your machine. In order to setting up the project, from the root of the project do:

```shell
mkdir build
cd build
cmake .. 
```

Then each application can be compiled by simply doing a

```shell
make <filename_without_extension>
```

If you run
```shell
make help
```
a list of all possible targets is displayed. 

Finally, as usual, running 

```shell
make clean
```

will delete all the executable and the compilation objects.


#### Code Documentation

In order to generate the documentation you need have doxygen and graphviz installed. On linux:
```shell
sudo apt-get install doxygen graphviz.
```

On Mac OSX with homebrew
```shell
brew install doxygen graphviz.
```

Then a 
```shell
make doc
```
will generate the documentation in the ``doc`` folder of your build.
        
