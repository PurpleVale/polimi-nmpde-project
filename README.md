### Organization of the code

The directory structure is as follows:

Project/</br>
├── bin/</br>
├── build/</br>
├── include/</br>
├── lib/</br>
├── src/</br>
└── test/</br>

- `bin` includes the binary output for all compilation
- `build` contains the output of CMAKE
- `incldue` contains all headers file
- `lib` contains the compiled libraries
- `src` contains all the C++ code to be compiled
- `test` contains only C++ code to be used for testing

### Compiling

If the code is executed in AppTainer run the followings:

```zsh
module load gcc-glibc dealii
```

To compile we use CMAKE, from the `Project` directory run:

```zsh
cmake -S ./ -B ./build/
cd build
make
cd ..
```

If you are using Clion as IDE then the .idea folder should automatically load options to compile if the libraries are installed on the machine.

### Executing

The Executable are located in `Project/bin` without any exension, simply execute them from the `Project` directory

```zsh
./bin/binary_name
```
