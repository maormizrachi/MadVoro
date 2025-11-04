# How to Cite This Project

If you use **MadVoro** in your research, publications, or derived software, please cite the following paper:

> **Mizrachi, M.**, **Raveh, B.**, & **Steinberg, E.** (2025).  
> *madvoro: parallel construction of Voronoi diagrams in distributed memory systems*.  
> **RAS Techniques and Instruments**, 4, rzaf039.  
> https://doi.org/10.1093/rasti/rzaf039

---

## BibTeX Citation

```bibtex
@article{10.1093/rasti/rzaf039,
    author  = {Mizrachi, Maor and Raveh, Barak and Steinberg, Elad},
    title   = {madvoro: parallel construction of Voronoi diagrams in distributed memory systems},
    journal = {RAS Techniques and Instruments},
    volume  = {4},
    pages   = {rzaf039},
    year    = {2025},
    month   = {09},
    issn    = {2752-8200},
    doi     = {10.1093/rasti/rzaf039},
    url     = {https://doi.org/10.1093/rasti/rzaf039},
    eprint  = {https://academic.oup.com/rasti/article-pdf/doi/10.1093/rasti/rzaf039/64231338/rzaf039.pdf},
    abstract = {Voronoi diagrams are essential geometrical structures with numerous applications, particularly astrophysics-driven finite volume methods. While serial algorithms for constructing these entities are well-established, parallel construction remains challenging. This is especially true in distributed memory systems, where each host manages only a subset of the input points. This process requires redistributing points across hosts and accurately computing the corresponding Voronoi cells. In this paper, we introduce a new distributed construction algorithm, which is implemented in our open-source C++ 3D Voronoi construction framework. Our approach leverages Delaunay triangulation as an intermediate step, which is then transformed into a Voronoi diagram. We introduce the algorithms we implemented for the precise construction and our load-balancing approach and compare the running time with other state-of-the-art frameworks. madvoro is a versatile tool that can be applied in various scientific domains, such as mesh decomposition, computational physics, chemistry, and machine learning.}
}
