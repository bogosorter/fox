import os
import subprocess
import json

print('compiling', end='')
subprocess.run(['make'], check = True)

times = dict()

for dirpath, _, filenames in os.walk('tests'):

    # Test small inputs first
    filenames.sort(key = lambda name: (len(name), name))

    for filename in filenames:
        if 'input' not in filename: continue
        
        n = int(filename[5:])

        input_path = f'{dirpath}/{filename}'
        with open(input_path) as f:
            input = f.read()

        print(f'\rtesting input {input_path}')

        times[n] = dict()
        for p in [1, 2, 3, 4, 5, 6]:
            if n % p != 0: continue

            times[n][p * p] = []
            for _ in range(2):
                result = subprocess.run(['mpirun', '-np', str(p * p), '--hostfile', 'hostfile.txt', './fox'], input = input, capture_output = True, text = True, check = True).stdout
                time = float(result.splitlines()[0])
                times[n][p * p].append(time)

with open('output.json', 'w') as f:
        json.dump(times, f, indent=4)
