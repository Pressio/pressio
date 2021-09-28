
# rom: Galerkin: General Info

\todo: write this better


At a high level, using the pressio Galerkin ROMs involves three main steps:

1. *create*: you create an instance of a "Galerkin problem"

2. *extract*: you extract the underlying stepper object owned by the problem

3. *solve*: you use the stepper to solve in time the Galerkin problem


You should now pause and think for a second about the steps above.
What does a stepper have to do with a Galerkin ROM?
The answer is that practically speaking, at the lowest-level,
a Galerkin problem can be reduced to simply a "custom" stepper to advance in time.
This is exactly how pressio implements this and the reason why a Galerkin
problem contains a stepper object inside: when you create the
problem, pressio creates the appropriate custom stepper
object that you can use. You don't need to know how this is done,
or rely on the details, because these are problem- and implementation-dependent,
and we reserve the right to change this in the future.


Pressio currently support three variants of a Galerkin problem:

- Default: [link](md_pages_components_rom_galerkin_default.html)
- Hyper-reduced: [link](md_pages_components_rom_galerkin_hypred.html)
- Masked: [link](md_pages_components_rom_galerkin_masked.html)
