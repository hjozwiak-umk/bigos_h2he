---
project: the SCATTERING code
summary: BIGOS – Binary Inelastic scattering and Generalized Optical cross Section package, vs. 0.01
src_dir: src
include: src
output_dir: doc
media_dir: media
display: public
         protected
source: true
proc_internals: true
md_extensions: markdown.extensions.toc
graph: true
graph_maxnodes: 300
graph_maxdepth: 5
coloured_edges: true
sort: permission-alpha
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html
print_creation_date: true
creation_date: %Y-%m-%d %H:%M %z
author: Hubert Jóźwiak
author_email: hubert.jozwiak@doktorant.umk.pl
github: https://github.com/hjozwiak-umk
project_github: https://github.com/hjozwiak-umk/bigos_h2he
project_download: https://github.com/hjozwiak-umk/bigos_h2he.git
dbg: true
page_dir: pages
---

**BIGOS** – **B**inary **I**nelastic scattering and **G**eneralized **O**ptical cross **S**ection package is a set of FORTRAN 90 codes that
allow the user to obtain various physical quantities related to the scattering problem involving diatomic molecules and
atoms. The package is being developed at the Nicolaus Copernicus University in Toruń.

![bigos_graph](|media|/bigos_graph.jpg "Structure of the BIGOS project")

Here, we present the SCATTERINC code, the central part of the BIGOS package.
The purpose of the SCATTERING code is to solve the coupled equations for a given scattering system, provide the
scattering S-matrix elements and calculate the state-to-state cross-sections.

**Please read the [Program Description](page/index.html).**

@note
This version of the code is adjusted for diatom - atom collision systems,<br>
in particular the H<sub>2</sub> - He system. <br>
**Please, refer to this version of the code by citing the following paper** <br>
H. Jozwiak, F. Thibault, A. Viel, P. Wcislo, F. Lique, <br>
Rovibrational (de-)excitation of H<sub>2</sub> by He revisited <br>
https://doi.org/10.48550/arXiv.2311.09890
