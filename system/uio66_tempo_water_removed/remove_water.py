from random import randint

count = 128   # count of water molecules
a_prev = 3535 # current count of atoms in .gro file (without water)

with open('water1041.gro', 'r') as f:
    lines = [line for line in f]
f.close()

mol_lines = []
for mol_idx in range(int(len(lines) / 4)):
    mol_lines.append([])
    for atom_idx in range(4):
        mol_lines[-1].append(lines[mol_idx * 4 + atom_idx])

selected = []
while len(selected) < count:
    val = randint(0, 1041)
    if val not in selected:
        selected.append(val)
selected.sort()

def format_val(val, cnt):
    val = str(val)
    while len(val) < cnt:
        val = ' ' + val
    return val

with open(f'water{count}.gro', 'w') as f:
    m_cnt = 0
    a_cnt = a_prev
    for idx in selected:
        m_cnt += 1
        for atom_idx in range(4):
            a_cnt += 1
            line = mol_lines[idx][atom_idx]
            line = format_val(m_cnt, 5) + line[5:]
            line = line[0:16] + format_val(a_cnt, 4) + line[20:]
            f.write(line)
f.close()



