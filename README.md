# Design and implementation of a Python module for quaternions
> Codes for Master thesis in Uni Stuttgart
---
![CoverPage](https://github.com/Goessi/Master_Thesis/blob/master/coverPage1.png)
### Table of Contents
- [Author Info] (#author-info)
- [Description] (#description)
- [How To Use] (#how-to-use)
- [License] (#license)
---
## Author Info
|  Author  |         Email          |
|----------|------------------------|
| Jing QIN | Jing.QIN94@outlook.com |
---
## Description
This module is developed for quaternions. By using this module, calculations between quaternions are achieved. Independent functions in functions.py are used for Direction-Cosine Matrix(DCM). 

#### Technologies and Modules
- Python 3.7
- Cython
- NumPy 1.16.2
- Decimal 1.70
---
## How To Use
> download Quaternion.py and functions.py to the same folder the test file belongs to.
#### Import
~~~~
from Quaternion import *
~~~~
Or
~~~~
import Quaternion
~~~~
#### Initialization
~~~~
q1 = Quaterinon(1,2,3,4)
q2 = Quaternion(1.1,1.2,1.3,1.4)
~~~~
#### Operations
1.Addition
~~~~
q1 + q2
~~~~
2.Subtraction
~~~~
q1 - q2
~~~~
3.Hamilton Product
~~~~
q1 * q2
~~~~
4.Dot product
~~~~
q1.dot()
~~~~
5.Conjugate
~~~~
q1.conj()
~~~~
6.Scale
~~~~
q1.scale_mul(scalar)
~~~~
7.Norm
~~~~
q1.norm()
~~~~
8.Normalization
~~~~
q1.norm_q()
~~~~
9.Space rotation
~~~~
q1.rotator(rotation angle in radians, axes)
~~~~
10. Transformation
~~~~
q1.toDCM()
~~~~
---
## License

---
