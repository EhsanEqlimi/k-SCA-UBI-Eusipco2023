# k-SCA-UBI-Eusipco2023

## 📄 Home

Welcome to the **UBI-(m−1)-Sparsity Solver**!

This repository provides the implementation for the method introduced in our paper:  
📘 *"Eqlimi E, Makkiabadi B, Kouti M, Fotouhi A, Sanei S. Underdetermined Blind Identification via k-Sparse Component Analysis: RANSAC-driven Orthogonal Subspace Search. in proceeding of EUSICO 2023, Helsinki, Finland"*.

Our algorithm addresses the challenge of **blind source separation** in the case where sources are sparse and the number of active sources at each time point equals *m–1*, where *m* is the number of sensors. This scenario often arises in real-world applications with high noise levels.

---

## 🚀 Overview

### Problem

**Underdetermined Blind Identification (UBI):**  
Recovering source signals when there are more sources than sensors.

### Traditional Approaches

- **SCA:** Assumes only 1 active source per time instant.  
- **k-SCA:** Allows up to *k* active sources but usually relies on 1-sparse assumptions.

### Challenge Addressed

Existing k-SCA methods struggle when *k = m−1*, especially under noisy conditions.

---

## 🧠 Our Contribution

We propose a novel two-step method that works efficiently when *k = m−1*:

1. **Estimate orthogonal complement subspaces** using a RANSAC-based sampling approach.
2. **Identify mixing vectors** using a Gram-Schmidt orthogonalization process.

This approach is computationally efficient and robust to noise.

---

Stay tuned for more detailed instructions, examples, and updates. If you use this work, please consider citing our paper.
