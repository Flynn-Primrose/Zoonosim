import Zoonosim as zn

print(zn.__path__)

zn_attr = dir(zn)

for attr in zn_attr:
    print(attr)