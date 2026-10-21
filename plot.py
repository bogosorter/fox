import matplotlib.pyplot as plt
import json
import math

with open('robin.json') as f:
    data = json.load(f)

def create_series(data, n):
    x = [int(key) for key in data[str(n)].keys()]
    y = [sum(row) / len(row) for row in data[str(n)].values()]
    return x, y


_, ax = plt.subplots()

x, y = create_series(data, 5)
plt.plot(x, y, label="5")

x, y = create_series(data, 6)
plt.plot(x, y, label="6")

x, y = create_series(data, 300)
plt.plot(x, y, label="300")

x, y = create_series(data, 900)
plt.plot(x, y, label="900")

x, y = create_series(data, 1200)
plt.plot(x, y, label="1200")

ax.legend()
plt.show()
