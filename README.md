# Shtawa
Numerical analysis dydactical library for modern Fortran. The goal of the library is to contain well-documented and well-tested functions and subroutines, with examples and notes about the theory behind each method and its limitations. The goal is not to produce the most efficient implementation.

## Dependencies
Library itself: gfortran compliant with the 2008+ standard. Some examples use Python (matplotlib, numpy) for plotting.

## Installation
Clone the repository or download it. Then, go inside the main directory (`cd path_to_shtawa_repo`) and compile the main library with `make`. This produces a shared object (.o) and a module file (.mod) in the directory `bin/`. To clear the compilation, do `make clean`.

Shtawa is composed of functions and subroutines ("procedures") which try to be dependency-free, so that their source code can be reused in other projects if needed. They are classified in directories by topic:

| Directory | Topic |
| ----------|-------------|
| `roots/`  | Root-finding algorithms (solutions to equations) |
| `linalg/` | Linear algebral (matrices, linear systems of equations) |

The in a Fortran source file, the procedures are accessible by a module
```use shtawa```
and including the binary objects during compilation (see the makefiles in the `examples/` directories for reference on how to link the library to your source code).


## Use and examples
The `examples/` directory contains examples of use of every function and subroutine in the library, organized by topic. To compile an example (e.g. `examples/roots/bisection_example.f90`) go to its directory (e.g. `examples/roots/`) and do `make example_name` (in this example, `make bisection_example`). This produces a binary called `example_name.bin` (e.g., `bisection_example.bin`) that can be run with `./example_name.bin`. See inside of each example for details on their inputs and outputs. To remove the binaries, do `make clean`. Some binaries create data, which is stored by default in a `data/` folder. To remove the data from all examples, do `make clean_data`

## Name
(Shtáwa) "To count" in Chirripo Cabécar language (http://www.ggc.cr/cabecar/diccionario/esp-cab/)

## License
BSD-0 License (essentially public domain)

Copyright (C) 2022 by André Oliva

Permission to use, copy, modify, and/or distribute this software for any purpose with or without fee is hereby granted.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.