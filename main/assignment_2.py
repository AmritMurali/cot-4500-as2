import numpy as np
np.set_printoptions(precision=7, suppress=True, linewidth=100)

# Assignment 2 for COT 4500
# Amrit Murali 

def neville(num, x, vals, w):
  nev = np.zeros((num, num))
  for i in range(num): # left-most column are values
    nev[i][0] = vals[i]

  for j in range(1, num):
    for k in range(1, j + 1):
      term1 = (w - x[j-k]) * nev[j][k-1]
      term2 = (w - x[j]) * nev[j-1][k-1]
      nev[j][k] = (term1-term2)/(x[j]-x[j-k])
  return nev[num-1][num-1]

def divdif(num, x, vals):
  dd = np.zeros((num, num))
  for i in range(num):
    dd[i][0] = vals[i]
  for j in range(1, num):
    for k in range(1, j + 1):
      dd[j][k] = (dd[j][k-1]-dd[j-1][k-1])/(x[j] - x[j-k])
  return dd


def hermite(num, x, f, fp):
  h = np.zeros((num * 2, num * 2))
  for i in range(num):
    h[i * 2][0] = x[i]
    h[i * 2 + 1][0] = x[i]
    h[i * 2][1] = f[i]
    h[i * 2 + 1][1] = f[i]
    h[i * 2 + 1][2] = fp[i]

  for j in range(2, num * 2):
    for k in range(j-1, num * 2):
      if(h[k][j] == 0 and h[k][j-1] != 0 and h[k-1][j-1] != 0):
        h[k][j] = (h[k][j-1]-h[k-1][j-1])/(h[k][0] - h[k-j+1][0])
  return h

def csi(n,x,y):
  A = np.zeros((n, n))
  h = np.zeros(n - 1)
  A[0][0] = 1
  A[n-1][n-1]=1
  for a in range(n-1):
    h[a] = x[a+1]-x[a]
  for b in range(0, n-2):
    A[b+1][b] = h[b]
    A[b+1][b+1] = 2*(h[b]+h[b+1])
    A[b+1][b+2] = h[b+1]
  print(A)
  c = np.zeros(n)
  for r in range(0, n-2):
    c[r+1]=(3/h[r+1])*(y[r+2]-y[r+1]) - (3/h[r])*(y[r+1]-y[r])
  print()
  print(c)
  print()
  z=np.dot(np.linalg.inv(A), c) 
  print(z)

# 1. using neville's method
# 2nd degree interpolating value for f(3.7)
# what is the y-value for this supposed x based on these data points?
num = 3 # data points
x = [3.6,3.8,3.9]
w = 3.7
vals = [1.675, 1.436, 1.318]
# if matrix is printed:
# [[1.675  0.     0.    ]
# [1.436  1.5555 0.    ]
# [1.318  1.554  1.555 ]]
# just the answer: 1.5549999999999995
print(neville(num, x, vals, w))
print()

# 2. using Newton's forward method
num = 4 # data points
x = [7.2,7.4,7.5,7.6]
vals = [23.5492, 25.3913, 26.8224, 27.4589]
# first find divided difference table
# delta y / delta x repeated
# helps find differences of higher order
dd = divdif(num, x, vals)
two = []
for u in range(1, num):
  two.append(dd[u][u])
print(two)
print()

# 3. polynomial approximations for degree 1,2,and 3 for f(7.3)
apr = 7.3
sum = dd[0][0]
for a in range(num - 1): 
  term = dd[a+1][a+1]
  for b in range(0, a + 1):
    term = term * (apr - x[b])
  sum = sum + term
print(sum)
print()

# 4. hermite interpolation
num = 3
x = [3.6,3.8,3.9]
f = [1.675, 1.436, 1.318]
fp = [-1.195, -1.188, -1.182]
print(hermite(num, x, f, fp))
print()

#5. cubic spline interpolation
n = 4
x = [2, 5, 8, 10]
y = [3, 5, 7, 9]
csi(n, x, y)
print()