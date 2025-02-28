def parse_rates(s):
    parts = s.split(" ")
    wr = []
    for part in parts[:-1]:
        part = part.strip("()")
        sub = part.split(",")
        wr.append((float(sub[0]), float(sub[1])))
    return wr

def E(wr):
    e = 0
    for weight, rate in wr:
        e += rate * weight
    return e

def var(wr):
    mu = E(wr)
    v = 0
    for weight, rate in wr:
        v += weight * (rate - mu) * (rate - mu)
    return v
