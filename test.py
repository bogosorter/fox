import os
import subprocess

print('compiling', end='')
subprocess.run(['gcc', '-Wall', 'fox.c', '-o', 'fox'], check = True)

failed = 0

for dirpath, _, filenames in os.walk('tests'):

    # Test small inputs first
    filenames.sort(key = lambda name: (len(name), name))

    for filename in filenames:
        if 'input' not in filename: continue

        input_path = f'{dirpath}/{filename}'
        expected_path = f'{dirpath}/{filename.replace('input', 'output')}'

        with open(input_path) as f:
            input = f.read()

        with open(expected_path) as f:
            expected = f.read()

        print(f'\rtesting input {input_path}', end='')
        result = subprocess.run('./fox', input = input, capture_output = True, text = True, check = True).stdout

        if result.strip() != expected.strip():
            failed += 1
            print(f'\rinput {input_path} failed')

if failed == 0: print('\r\033[Kpassed all tests')
else: print(f'\r\033[Kfailed {failed} tests')
