project: quadpack
src_dir: ./src
output_dir: ./doc
media_dir: ./media
project_github: https://github.com/jacobwilliams/quadpack
summary: Modern Fortran QUADPACK
author: Jacob Williams
github: https://github.com/jacobwilliams
predocmark_alt: >
predocmark: <
docmark_alt:
docmark: !
display: public
source: true
graph: true
exclude: **/quadpack_single.F90
         **/quadpack_double.F90
         **/quadpack_quad.F90
         **/quadpack.F90
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html

{!README.md!}