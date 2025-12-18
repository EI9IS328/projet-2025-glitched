import glob
import re
import struct
import os

FILEMAGIC = 0xCAFEBABE
FILEMAGIC_FORMAT = "<I"
FILEMAGIC_SIZE = 4
HEADER_FORMAT = "<iiii"
HEADER_SIZE = struct.calcsize(HEADER_FORMAT)
ROW_FORMAT = "<ffff"
ROW_SIZE = struct.calcsize(ROW_FORMAT)

def openBin():
    file = "snapshot600.bin"
    dir = os.path.dirname(__file__)
    path = os.path.normpath(os.path.join(dir, '..', '..', "snapshot", file))

    with open(path, "rb") as f:
        sig = struct.unpack(HEADER_FORMAT, f.read(HEADER_SIZE))
        print(sig)


if __name__ == "__main__":
    openBin()