import random
import math
import matplotlib.pyplot as plt

N = 1000
N_inside = 0

x_values_inside = []
y_values_inside = []

x_values_outside = []
y_values_outside = []

with open("points_outside.txt","w") as file_outside, open("points_inside.txt","w") as file_inside, open("convergence.txt", "w") as file_convergence, open("error.txt","w") as file_error:
    for i in range(1,N):
        x = random.random()
        y = random.random()
        if(x*x+y*y <= 1):
            N_inside = N_inside + 1
            file_inside.write(f"{x} {y}\n")
            x_values_inside.append(x)
            y_values_inside.append(y)
        else:
            file_outside.write(f"{x} {y}\n")
            x_values_outside.append(x)
            y_values_outside.append(y)
        pi = 4*N_inside/i
        file_convergence.write(f"{i} {pi}\n")
        error = abs(pi/math.pi - 1)
        file_error.write(f"{i} {error}\n")

print("The irrational number pi is: ", pi)

plt.scatter(x_values_inside, y_values_inside, color='blue', marker='.')
plt.scatter(x_values_outside, y_values_outside, color='green', marker='.' )
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Randomly Generated points')
plt.grid(True)
plt.show()

