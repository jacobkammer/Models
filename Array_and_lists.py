import numpy as np

for i in range(0,10):
     lista = []
     lista.append(i)
print(lista)
print(type(lista))

listb = []
for i in range(0,10):
    listb.append(i)
print(listb)
print(type(listb))



numbers = np.arange(9)
numbers_2d = numbers.reshape(3, 3)
#add a row to the array.  single arguement in the form of a tuple
new_row = np.array([34, 78, 12])
results02 = np.vstack((numbers_2d, new_row))
print(results02)
print(type(results02))



A = np.array([[1, 2, 3], [45, 4, 7], [1, 3, 89]])
#add a column to the array. single arguement in the form of a tuple
new_column = np.array([56, 46, 78])
result = np.column_stack((A, new_column))
print(result)
print(type(result))

trig = []
B = [0, 3*np.pi/4, np.pi/2, 2*np.pi, 2*np.pi/3, np.pi, -np.pi, -np.pi/2]
for i in B:
    a = np.cos(i)
    trig.append(a)
h = zip(B, trig)
print(dict(h))
print(type(dict(h)))

print(type(h))


#notes on zip function. combines multiple iterables into one set of elements a zip object, zip object is a tuple
#list comprehension.  Equivalent to the above for loop
cosines = [np.cos(i) for i in B]
print(cosines)


#another list comphrehension
#pattern:  (values) = [(expression) for (value) in (collection)]
squares = []
for x in range(10):
    squares.append(x*x)
print(squares)

squares = [x*x for x in range(10)]
print(squares)

import numpy as np
array = np.arange(10)
print(array)