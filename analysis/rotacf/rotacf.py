import math

def read_xvg(file_name):
    file = open(file_name, 'r')
    lines = [line.replace('\n', '') for line in file]
    lines = [line for line in lines if len(line) != 0 \
             if line[0] != '@' and line[0] != '#' and line[0] != '&']
    file.close()

    x = [float(line.split()[0]) for line in lines]
    y = [float(line.split()[1]) for line in lines]
    return x, y

def calculate_tau(x, y):
    for val_idx in range(len(x)):
        if y[val_idx] < y[0] / 3.14:
            return x[val_idx] - x[0]
    return (x[-1] - x[0]) / math.log(y[0] / y[-1])

x, y = read_xvg('/home/sc_enigma/Projects/MD_UIO66/analysis/rotacf/N_O/uio66_tempo.xvg')
print(f'Tau = {calculate_tau(x, y) / 1000.0} ns')