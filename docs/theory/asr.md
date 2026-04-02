# **Acoustic Sum Rules and Invariance Conditions**

The force constant matrix $\Phi_{i \alpha j \beta}(\mathbf{R})$ describes the interaction between atom $i$ and atom $j$ in a unit cell displaced by $\mathbf{R}$. It is defined as the second derivative of the potential energy with respect to atomic displacements.
$$
\Phi_{i \alpha j \beta}(\mathbf{R}) = \Phi_{j \beta i \alpha }(-\mathbf{R}) 
$$

## 1. **Acoustic Sum Rules**
Below are Acoustic Sum Rules (ASR) obeyed by the force constant matrix. These rules are derived from the invariance of the total energy under global translations and rotations.

### **Translation Invariance**:
This condition ensures that a rigid translation of the entire crystal does not change its potential energy.
$$
\sum_{\mathbf{R}, j} \Phi_{i \alpha j \beta}(\mathbf{R}) = 0,
$$

where $i,j$ denote the atomic indices and $\alpha,\beta,\gamma,\delta$ denote the polarization directions. $\mathbf{R}$ is the position vector of the unit cell. The summation over $\mathbf{R}$ is in the Wigner-Seitz cell. 

### **Rotation Invariance**:
This condition ensures that a rigid rotation of the entire crystal does not change its potential energy.
$$
\sum_{\mathbf{R}, j} \left( \Phi_{ i \alpha j \beta} (\mathbf{R}) \tau^{\gamma}_{\mathbf{R} j} - \Phi_{ i \alpha j \gamma} (\mathbf{R}) \tau^{\beta}_{\mathbf{R} j} \right) = 0,
$$

where $\boldsymbol{\tau}_{\mathbf{R} j} = \mathbf{R} + \boldsymbol{\tau}_j$ is the position of atom $j$ in the unit cell located at $\mathbf{R}$. $\boldsymbol{\tau}_j$ is the position of atom $j$ with respect to the unit cell.


### **Huang Conditions**:
Huang conditions are related to the symmetry of the elastic constant tensor and are necessary for the consistency between the microscopic force constants and macroscopic elasticity.
$$
\sum_{\mathbf{R}, i, j} \Phi_{i \alpha j \beta} (\mathbf{R}) \tau^{\gamma}_{\mathbf{R} i j} \tau^{\delta}_{\mathbf{R} i j} = \sum_{\mathbf{R}, i, j} \Phi_{ i \gamma j \delta} (\mathbf{R}) \tau^{\alpha}_{\mathbf{R} ij} \tau^{\beta}_{\mathbf{R} ij}
$$

where $\boldsymbol{\tau}_{\mathbf{R} i j} = \boldsymbol{\tau}_{i} - \boldsymbol{\tau}_{j} - \mathbf{R}$ is the vector between two atoms.

---

## 2. **Enforcing Acoustic Sum Rules**

It is often the case that the force constants obtained from standard *ab initio* calculations break acoustic sum rules due to the finite size of the grid or basis set, which leads to non-zero frequencies for the acoustic modes at the $\Gamma$ point.

In order to enforce ASR, the first step is to explicitly write the invariance conditions in the form of a matrix-vector product $A\Phi = 0$, where $A$ is the constraint matrix and $\Phi$ is the vectorized force constant matrix. Here, $N_A$ refers to the number of atoms in the unit cell.

### 1. Translation Operator

The translation operator tensor `A` is defined as:

$$
A^{i\alpha\beta}_{\mathbf{R}i'\alpha'j'\beta'} = \delta_{ii'} \delta_{\alpha\alpha'} \delta_{\beta\beta'} 
$$

**Constraint:**

$$
\sum_{\mathbf{R}, i', \alpha', j', \beta'} A^{i \alpha \beta}_{\mathbf{R} i' \alpha' j' \beta'} \Phi_{\mathbf{R} i' \alpha' j' \beta'} = 0 \quad \Rightarrow \quad 9 N_A
$$

### 2. Rotation Operator

The rotation operator tensor `A` is:

$$
A^{i\alpha\beta\gamma}_{\mathbf{R}i'\alpha'j'\beta'} = \delta_{ii'} \delta_{\alpha\alpha'} \left( \delta_{\beta\beta'} \tau^{\gamma}_{\mathbf{R}j'} - \delta_{\gamma\beta'} \tau^{\beta}_{\mathbf{R}j'} \right)
$$

**Symmetry & Constraints:**

$$
A^{i\alpha\beta\gamma} = - A^{i\alpha\gamma\beta} \quad \longrightarrow \quad 3 \cdot N_A \cdot (3) = 9 N_A
$$

**Resulting Equation ($A\Phi=0$)**:

$$
\sum_{\mathbf{R}, j'} \left( \Phi_{\mathbf{R} i \alpha j' \beta} \tau^{\gamma}_{\mathbf{R} j'} - \Phi_{\mathbf{R} i \alpha j' \gamma} \tau^{\beta}_{\mathbf{R} j'} \right) = 0
$$



### 3. Huang Operator

The Huang operator tensor `A` is:

$$
A^{\alpha\beta\gamma\delta}_{\mathbf{R}i'\alpha'j'\beta'} = \left[ \delta_{\alpha\alpha'} \delta_{\beta\beta'} \tau^{\gamma}_{\mathbf{R}i'j'} \tau^{\delta}_{\mathbf{R}i'j'} - \tau^{\alpha}_{\mathbf{R}i'j'} \tau^{\beta}_{\mathbf{R}i'j'} \delta_{\gamma\alpha'} \delta_{\delta\beta'} \right]
$$

**Symmetry & Constraints:**

For $(\alpha, \beta) \longleftrightarrow (\gamma, \delta)$, $A^{\alpha\beta\gamma\delta}$ is anti-symmetric i.e,

$$
A^{\alpha\beta\gamma\delta} = - A^{\gamma\delta\alpha\beta}
$$


For $(\alpha \leftrightarrow \beta)$ and $(\gamma \leftrightarrow \delta)$ exchanges, $A^{\alpha\beta\gamma\delta}$ is symmetric.

$$
\left. \begin{matrix} (\alpha \leftrightarrow \beta) \\ (\gamma \leftrightarrow \delta) \end{matrix} \right\} \longrightarrow \underline{\underline{15}}
$$

Total constraints without reduction : 81 
* Due to anti-symmetry property for $(\alpha, \beta) \leftrightarrow (\gamma, \delta)$, it reduces to 36 constraints.
* Due to symmetry property for $(\alpha \leftrightarrow \beta)$ and $(\gamma \leftrightarrow \delta)$ further reduces to 15 constants.

---



