f = open("runs.txt", "r");
lines = 0;
max = 0;
i = 0;
for line in f:
    try:
        if (float(line.split()[6]) > max and float(line.split()[6]) < 100):
            max = float(line.split()[6]);
            lines = i;
    except:
        None
    i += 1;
print(lines);
print(max);