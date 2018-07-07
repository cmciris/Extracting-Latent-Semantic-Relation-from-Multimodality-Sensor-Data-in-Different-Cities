=================
I. Contents 
=================
1. Non-convex regularized bundle method (NRBM version 0.7b) for machine learning.
2. OCR data in Matlab format
3. Implementation of linear chain conditional random fields (CRFs) with different training criteria : 
	- conditional likelihood (Lafferty et al. 2001)
	- max-margin (Taskar et al. 2004)
	- PosLearn (Sarawagi et al. 2008)
4. Implementation of linear chain neural conditional random fields (NeuroCRFs) (Do et Artieres, AISTATS 2010)

=================
II. Install
=================
Under Matlab, launch makeNRBM at NRBM root folder

=================
III. Usage
=================
See example_crf.m and example_neurocrf.m
Ex. Under matlab, launch exmaple_crf at NRBM root folder

=================
IV. Requirements
=================
SVM-like quadratic programming solvers
1. libqp (http://cmp.felk.cvut.cz/~xfrancv/libqp/html/) 
2. stprtool (http://cmp.felk.cvut.cz/cmp/software/stprtool/). 
These libraries are included in this distribution of NRBM for usage simplicity, but their availability are NOT guaranteed and subject to removal. 

=================
V. References
=================

[1] Trinh-Minh-Tri Do and Thierry Artieres. Regularized Bundle Methods for Convex and Non-Convex Risks. Accepted for publication, JMLR 2013.

[2] Trinh-Minh-Tri Do and Thierry Artieres. Large Margin Training for Hidden Markov Models with Partially Observed States. ICML 2009.

[3] Trinh-Minh-Tri Do and Thierry Artieres. Neural conditional random fields. AISTATS 2010.

[4] Trinh-Minh-Tri Do. Regularized bundle methods for large-scale learning problems with an application to large margin training of hidden Markov models. Phd dissertation, Pierre and Marie Curie University, Paris, 2010.


=================
VI. License
=================

Copyright (C) 2012, by Trinh-Minh-Tri Do minhtrido@gmail.com

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


=== last revision: Nov-23-2012 ===