# The Schwarz-Christoffel transformation

A conformal map $f: \mathbb{C} \rightarrow \mathbb{C}$ is an angle-preserving transformation of points located in the complex plane.
Conformal maps are particularly useful for solving electrostatics problems governed by Laplace's equation since solutions $w$ are harmonic and conformal maps preserve this property.
Hence, if an electrostatics boundary value problem can be solved in one geometry, it can be mapped onto a different boundary geometry.
The Schwarz-Christoffel transformation implements the mapping from a circular boundary to an arbitrary polygon.

This library implements

$$w = f(z) = c \int_0^z \prod_{k=1}^N \left(1-\frac{\tilde{z}}{z_k}\right)^{-\beta_k}\,{\rm d}\tilde{z} \ , $$

where $c$ is a constant for scaling and rotating, and $z_k \in \mathbb{C}$ are so-called prevertices which map to the vertices $w_k = f(z_k)$ of the polygon.
Finding the prevertices to a given polygon is called the Schwarz-Christoffel parameter problem and is also implemented.

<img width="163" alt="{21CA075D-98A8-4ABF-BC9E-3E83152B65DC}" src="https://github.com/user-attachments/assets/e0e1446c-b104-452e-a2ec-544c49cb4874" />
<img width="161" alt="{CEABC351-E375-4B62-A6E5-B088AE5D8EC4}" src="https://github.com/user-attachments/assets/b0fbf867-4529-4d7b-968e-7f136c346340" />

Vertices at infinity are also supported, enabling open boundary geometries.

<img width="161" alt="{762B762F-8609-4E7F-B5FA-8FBA4166926B}" src="https://github.com/user-attachments/assets/3566278b-4e9b-4ff9-be82-894bdc5acdea" />
<img width="163" alt="{B41AAB40-0B7B-4499-A643-4A72014EF8DC}" src="https://github.com/user-attachments/assets/4cc831e9-f395-4818-bebe-6ad338d02882" />
