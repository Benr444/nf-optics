# DOE-Design Program.

Author: Benjamin Ranson (bengr444@gmail.com)

---

## Setup.

###  1: Install *Julia*.

[*Julia*](https://julialang.org/) is a Python-like general-purposes programming language that's increasingly popular in scientific-computing applications due to its intuitive, mathlike syntax and high computation speeds, which are the reasons I chose it. I think you will find *Julia* enjoyable and perhaps you will consider using it on your own projects.

Go to https://julialang.org/downloads/ and follow the installation instructions there. Make sure *Julia* is added to your [PATH environmental variable](https://john-dugan.com/path-variable-in-windows/).

### 2: Test *Julia*.

Open your terminal of choice (on Linux) or `cmd` on windows.

Type `julia` and press enter. The *Julia* prompt should appear. *Julia* can be used as a line-by-line terminal execution (called the REPL) or it can execute script files.

Try running `x = (1:5).^2` in the terminal. You should get

	5-element Vector{Int64}:
	1
	4
	9
	16
	25

which is a list of the square of all the numbers 1 to 5. *Julia*'s convenient '.'-syntax to 'broadcast' math operations over collections is very useful for keeping code fast and concise.

To run a *Julia* script, write any amount of *Julia* code in a text file (conventionally with the file extension `.jl`.), then write `include("my-file.jl")` in your REPL.

Call `exit()` to stop the REPL and return to your main terminal.

### 3: Run DesignKinoform.jl.

To run the principle program here, simply write `include(DesignKinoform.jl)` in your REPL.

### 4: Adding missing packages.

You may get missing packages/libraries errors, however, the terminal/REPL should tell you how to fix them.

Another way to fix mising packages is, while at the REPL, type `]` (right square-bracket). Your prompt should change from `julia>` to `(@v1.X) pkg>` where `1.X` is your *Julia* version.

Once you have entered 'package mode' this way, you can simply type `add PackageName` and press enter to install a package and all its dependencies.

Typing `add Revise Images FileIO FFTW Random LaTeXStrings` shuld add all packages used in this repository at the time of writing.

You can then use 'control + c' to return to the REPL.

---

## Code Structure and Output.

The code is split among `DesignKinoform.jl` and `ROptics.jl`. Utility functions are put in the later, while problem-specific procedures are in the former.

The program will place its output in `out/`, which includes various graphs and raw matrix images. You can uncomment additional graphs in the source code to get more information, but it will slow down the runtime.

I've tried to comment the code well, but ultimately you should read `explanatory-paper.pdf` to understand its function and output.

## Using an IDE.

For a faster workflow, install an IDE like [VSCode](https://code.visualstudio.com/) its respective *Julia* extension. This can help manage the different files and terminals and aid in your editing with features like syntax-highlighting. Eventually, you will want to do this if you plan on editing the source.