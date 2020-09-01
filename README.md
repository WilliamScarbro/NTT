# Optimizing NTT using Spiral for applications in Post Quantum Cryptography

## Internship Overview
* Gives mathematic foundation for NTT
* Describes our approach in adapting Spiral for NTT
* Incorperates one NTT optimization ( \psi powers absorbtion)

## Spiral Code
* ntt.g: first NTT implementation in Spiral
* ntt_nw.g: includes \psi powers for negative wrapped convolution
* See [Spiral User's Manual](https://www.spiral.net/doc/usermanual/) to run

## Sage Code
* MatSPL.ipynb: Simple demonstration of NTT in Spiral's MatSPL framework
* Run at [https://cocalc.com/](https://cocalc.com/)

## Resources
* For an in depth look at [Spiral](https://ieeexplore.ieee.org/document/8510983)
* NTT optimizations [Roy et al.](https://www.iacr.org/archive/ches2014/87310183/87310183.pdf), [Poppelmann et al.](https://drive.google.com/drive/u/0/folders/1fyT38SmWpSqVXkr2eXH-haAuA-NhMzLM)
