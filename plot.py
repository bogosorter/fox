import matplotlib.pyplot as plt
import json
import math

with open('data/all.json') as f:
    data = json.load(f)

with open('data/robin.json') as f:
    robin = json.load(f)

def create_series(data, n):
    x = [int(key) for key in data[str(n)].keys()]
    y = [sum(row) / len(row) for row in data[str(n)].values()]
    return x, y

def create_speedup(data, n):
    x = [int(key) for key in data[str(n)].keys()]
    y = [sum(row) / len(row) for row in data[str(n)].values()]
    base = y[0]
    y = [base / v for v in y]
    return x, y


_, ax = plt.subplots()

x, y = create_series(data, 5)
plt.plot(x, y, label = "5", marker = "o")

x, y = create_series(data, 6)
plt.plot(x, y, label = "6", marker = "o")

x, y = create_series(data, 300)
plt.plot(x, y, label = "300", marker = "o")

x, y = create_series(data, 900)
plt.plot(x, y, label = "900", marker = "o")

x, y = create_series(data, 1200)
plt.plot(x, y, label = "1200", marker = "o")

ax.set_ylim(bottom=0)
ax.legend()
plt.savefig("plots/time.png")

plt.close()

x, y = create_speedup(data, 5)
plt.plot(x, y, label = "5", marker = "o")

x, y = create_speedup(data, 6)
plt.plot(x, y, label = "6", marker = "o")

x, y = create_speedup(data, 300)
plt.plot(x, y, label = "300", marker = "o")

x, y = create_speedup(data, 900)
plt.plot(x, y, label = "900", marker = "o")

x, y = create_speedup(data, 1200)
plt.plot(x, y, label = "1200", marker = "o")

plt.savefig("plots/speedup.png")

plt.close()

x, y = create_series(robin, 5)
plt.plot(x, y, label = "5", marker = "o")

x, y = create_series(robin, 6)
plt.plot(x, y, label = "6", marker = "o")

x, y = create_series(robin, 300)
plt.plot(x, y, label = "300", marker = "o")

x, y = create_series(robin, 900)
plt.plot(x, y, label = "900", marker = "o")

x, y = create_series(robin, 1200)
plt.plot(x, y, label = "1200", marker = "o")

ax.set_ylim(bottom=0)
ax.legend()
plt.savefig("plots/time-robin.png")

plt.close()

x, y = create_speedup(robin, 5)
plt.plot(x, y, label = "5", marker = "o")

x, y = create_speedup(robin, 6)
plt.plot(x, y, label = "6", marker = "o")

x, y = create_speedup(robin, 300)
plt.plot(x, y, label = "300", marker = "o")

x, y = create_speedup(robin, 900)
plt.plot(x, y, label = "900", marker = "o")

x, y = create_speedup(robin, 1200)
plt.plot(x, y, label = "1200", marker = "o")

plt.savefig("plots/speedup-robin.png")
